// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <parallel/algorithm>
#include <omp.h>
#include <cassert>

#include "parallel.h"
#include "blockRadixSort.h"
#include "quickSort.h"
#include "utils.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"
#include "sliding_queue.h"
#include "dbg.h"

using namespace std;

typedef pair<uintE,uintE> intPair;
typedef pair<uintE, pair<uintE,intE> > intTriple;

template <class E>
struct pairFirstCmp {
  bool operator() (pair<uintE,E> a, pair<uintE,E> b) {
    return a.first < b.first; }
};

template <class E>
struct getFirst {uintE operator() (pair<uintE,E> a) {return a.first;} };

template <class IntType>
struct pairBothCmp {
  bool operator() (pair<uintE,IntType> a, pair<uintE,IntType> b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  }
};

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long n; // total number of characters
  char* Chars;  // array storing all strings
  long m; // number of substrings
  char** Strings; // pointers to strings (all should be null terminated)
  words() {}
words(char* C, long nn, char** S, long mm)
: Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {free(Chars); free(Strings);}
};

inline bool isSpace(char c) {
  switch (c)  {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ' : return true;
  default : return false;
  }
}

_seq<char> mmapStringFromFile(const char *filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG (sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char *p = static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
//  char *bytes = newA(char, n);
//  parallel_for(size_t i=0; i<n; i++) {
//    bytes[i] = p[i];
//  }
//  if (munmap(p, sb.st_size) == -1) {
//    perror("munmap");
//    exit(-1);
//  }
//  cout << "mmapped" << endl;
//  free(bytes);
//  exit(0);
  return _seq<char>(p, n);
}

_seq<char> readStringFromFile(char *fileName) {
  ifstream file (fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg (0, ios::beg);
  long n = end - file.tellg();
  char* bytes = newA(char,n+1);
  assert(bytes != NULL && "Malloc failure\n");
  file.read (bytes,n);
  file.close();
  return _seq<char>(bytes,n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, long n) {
  {parallel_for (long i=0; i < n; i++)
      if (isSpace(Str[i])) Str[i] = 0; }

  // mark start of words
  bool *FL = newA(bool,n);
  assert(FL != NULL && "Malloc failure\n");
  FL[0] = Str[0];
  {parallel_for (long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];}

  // offset for each start of word
  _seq<long> Off = sequence::packIndex<long>(FL, n);
  free(FL);
  long m = Off.n;
  long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char*, m);
  assert(SA != NULL && "Malloc failure\n");
  {parallel_for (long j=0; j < m; j++) SA[j] = Str+offsets[j];}

  free(offsets); 
  return words(Str,n,SA,m);
}

template <class vertex>
graph<vertex> readGraphFromFile(char* fname, bool isSymmetric, bool mmap) {
#ifndef VGR
  words W;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    assert(bytes != NULL && "Malloc failure\n");
    // Cannot mutate the graph unless we copy.
    parallel_for(size_t i=0; i<S.n; i++) {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    _seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyGraph") {
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph") {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
#ifndef WEIGHTED
  if (len != n + m + 2) {
#else
  if (len != n + 2*m + 2) {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  uintT* offsets = newA(uintT,n);
  assert(offsets != NULL && "Malloc failure\n");
#ifndef WEIGHTED
  uintE* edges = newA(uintE,m);
#else
  intE* edges = newA(intE,2*m);
#endif
  assert(edges != NULL && "Malloc failure\n");

  {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long i=0; i<m; i++) {
#ifndef WEIGHTED
      edges[i] = atol(W.Strings[i+n+3]);
#else
      edges[2*i] = atol(W.Strings[i+n+3]);
      edges[2*i+1] = atol(W.Strings[i+n+m+3]);
#endif
    }}
    
  //W.del(); // to deal with performance bug in malloc
  W.del(); //The original code ^ commented this out
#else // #ifdef VGR // Added by Priyank
  if ( sizeof(uintT) != sizeof(unsigned long long) ) {
      std::cout << sizeof(uintT) << " " << sizeof(unsigned long long) << std::endl;
      abort();
  }
  if ( sizeof(uintE) != sizeof(unsigned int) ) {
      std::cout << sizeof(uintE) << " " << sizeof(unsigned int) << std::endl;
      abort();
  }
  ifstream ifs(fname, std::ios::binary);
  if ( !ifs.good() ) {
	  std::cout << "Unable to open file: " << fname << std::endl;
	  abort();
  }
  unsigned long long major_version;
  unsigned long long minor_version;
  unsigned long long n1;
  unsigned long long m1;
  ifs.read((char*)&major_version, sizeof(unsigned long long));
  ifs.read((char*)&minor_version, sizeof(unsigned long long));
  ifs.read((char*)&n1, sizeof(unsigned long long));
  ifs.read((char*)&m1, sizeof(unsigned long long));

  long n,m;
  m = m1;
  n = n1;
#ifndef WEIGHTED
  if ( ( major_version != 1 ) || ( minor_version != 0 ) ) {
      std::cout << "major: " << major_version << " minor: " << minor_version << "  n: " << n << " m: " << m << std::endl;
	  abort();
  }
#else
  if ( ( major_version != 1 ) || ( minor_version != 4 ) ) {
      std::cout << "major: " << major_version << " minor: " << minor_version << "  n: " << n << " m: " << m << std::endl;
	  abort();
  }
#endif

	uintT* offsets = newA(uintT,n);
	assert(offsets != NULL && "Malloc failure\n");
	ifs.read((char*)&offsets[1], sizeof(unsigned long long) * (n-1));
	if ( !ifs || (ifs.gcount() != (sizeof(unsigned long long) * (n-1))) ) {
		std::cout << "Error in reading offsets." << std::endl;
		abort();
	}
	unsigned long long temp;
	ifs.read((char*)&temp, sizeof(unsigned long long) * (1));
	offsets[0] = 0;

#ifndef WEIGHTED
	uintE* edges = newA(uintE,m);
#else
	intE* edges_ = newA(intE,2*m);
	assert(edges_ != NULL && "Malloc failure\n");
	intE* edges = newA(intE,2*m);
#endif
	assert(edges != NULL && "Malloc failure\n");
#ifndef WEIGHTED
	ifs.read((char*)edges, sizeof(unsigned int) * m);
	if ( !ifs || (ifs.gcount() != (sizeof(unsigned int) * m)) ) {
		std::cout << "Error in reading edges." << std::endl;
		abort();
	}
#else
	ifs.read((char*)&edges_[0], sizeof(unsigned int) * m);
	if ( !ifs || (ifs.gcount() != (sizeof(unsigned int) * m)) ) {
		std::cout << "Error in reading edges." << std::endl;
		abort();
	}
    if ( (m % 2) == 1 ) {
        unsigned int x = 0;
        ifs.read((char*)&x, sizeof(unsigned int));
        //std::cout << "extra element: " << x << std::endl;
        assert(x == 0);
    }
	ifs.read((char*)&edges_[m], sizeof(unsigned int) * m);
	if ( !ifs || (ifs.gcount() != (sizeof(unsigned int) * m)) ) {
		std::cout << "Error in reading edges." << std::endl;
		abort();
	}
    //std::cout << edges_[m] << std::endl;
    {parallel_for(long i=0; i<m; i++) {
      edges[2*i] = edges_[i];
      edges[2*i+1] = edges_[i+m];
    }}
    free(edges_);
#endif

	ifs.close();
#endif // VGR

  vertex* v = newA(vertex,n);
  assert(v != NULL && "Malloc failure\n");

  {parallel_for (uintT i=0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l);
#ifndef WEIGHTED
    v[i].setOutNeighbors(edges+o);
#else
    v[i].setOutNeighbors(edges+2*o);
#endif
  }}

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    assert(tOffsets != NULL && "Malloc failure\n");
    {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair,m);
#else
    intTriple* temp = newA(intTriple,m);
#endif
    assert(temp != NULL && "Malloc failure\n");
    {parallel_for(long i=0;i<n;i++){
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsets);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<uintE>());
#else
    quickSort(temp,m,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<intPair>());
#else
    quickSort(temp,m,pairFirstCmp<intPair>());
#endif
#endif

    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE,2*m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    assert(inEdges != NULL && "Malloc failure\n");
    {parallel_for(long i=1;i<m;i++) {
#ifndef WEIGHTED
      inEdges[i] = temp[i].second;
#else
      inEdges[2*i] = temp[i].second.first;
      inEdges[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);

    {parallel_for(long i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
#ifndef WEIGHTED
      v[i].setInNeighbors(inEdges+o);
#else
      v[i].setInNeighbors(inEdges+2*o);
#endif
      }}

    free(tOffsets);
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges,inEdges);
    std::cout << "Read directed graph. Num Nodes = " << n << " and Num Edges = " << m << "\n";
    return graph<vertex>(v,n,m,mem);
  }
  else {
    free(offsets);
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges);
    std::cout << "Read undirected graph. Num Nodes = " << n << " and Num Edges = " << m << "\n";
    return graph<vertex>(v,n,m,mem);
  }
}

template <class vertex>
graph<vertex> readGraphFromBinary(char* iFile, bool isSymmetric) {
  char* config = (char*) ".config";
  char* adj = (char*) ".adj";
  char* idx = (char*) ".idx";
  char configFile[strlen(iFile)+strlen(config)+1];
  char adjFile[strlen(iFile)+strlen(adj)+1];
  char idxFile[strlen(iFile)+strlen(idx)+1];
  *configFile = *adjFile = *idxFile = '\0';
  strcat(configFile,iFile);
  strcat(adjFile,iFile);
  strcat(idxFile,iFile);
  strcat(configFile,config);
  strcat(adjFile,adj);
  strcat(idxFile,idx);

  ifstream in(configFile, ifstream::in);
  long n;
  in >> n;
  in.close();

  ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints
  in2.seekg(0, ios::end);
  long size = in2.tellg();
  in2.seekg(0);
#ifdef WEIGHTED
  long m = size/(2*sizeof(uint));
#else
  long m = size/sizeof(uint);
#endif
  char* s = (char *) malloc(size);
  in2.read(s,size);
  in2.close();
  uintE* edges = (uintE*) s;

  ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if(n != size/sizeof(intT)) { cout << "File size wrong\n"; abort(); }

  char* t = (char *) malloc(size);
  in3.read(t,size);
  in3.close();
  uintT* offsets = (uintT*) t;

  vertex* v = newA(vertex,n);
#ifdef WEIGHTED
  intE* edgesAndWeights = newA(intE,2*m);
  {parallel_for(long i=0;i<m;i++) {
    edgesAndWeights[2*i] = edges[i];
    edgesAndWeights[2*i+1] = edges[i+m];
    }}
  //free(edges);
#endif
  {parallel_for(long i=0;i<n;i++) {
    uintT o = offsets[i];
    uintT l = ((i==n-1) ? m : offsets[i+1])-offsets[i];
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      v[i].setOutNeighbors((uintE*)edges+o);
#else
      v[i].setOutNeighbors(edgesAndWeights+2*o);
#endif
  }}

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair,m);
#else
    intTriple* temp = newA(intTriple,m);
#endif
    {parallel_for(intT i=0;i<n;i++){
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
#ifndef WEIGHTED
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
#else
	temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
#endif
      }
      }}
    free(offsets);
#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<uintE>());
#else
    quickSort(temp,m,pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp,m,n+1,getFirst<intPair>());
#else
    quickSort(temp,m,pairFirstCmp<intPair>());
#endif
#endif
    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE,2*m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {parallel_for(long i=1;i<m;i++) {
#ifndef WEIGHTED
      inEdges[i] = temp[i].second;
#else
      inEdges[2*i] = temp[i].second.first;
      inEdges[2*i+1] = temp[i].second.second;
#endif
      if(temp[i].first != temp[i-1].first) {
	tOffsets[temp[i].first] = i;
      }
      }}
    free(temp);
    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);
    {parallel_for(long i=0;i<n;i++){
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
#ifndef WEIGHTED
      v[i].setInNeighbors((uintE*)inEdges+o);
#else
      v[i].setInNeighbors((intE*)(inEdges+2*o));
#endif
      }}
    free(tOffsets);
#ifndef WEIGHTED
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges,inEdges);
    return graph<vertex>(v,n,m,mem);
#else
    Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights,inEdges);
    return graph<vertex>(v,n,m,mem);
#endif
  }
  free(offsets);
#ifndef WEIGHTED
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges);
  return graph<vertex>(v,n,m,mem);
#else
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edgesAndWeights);
  return graph<vertex>(v,n,m,mem);
#endif
}

template <class vertex>
graph<vertex> readGraph(char* iFile, bool compressed, bool symmetric, bool binary, bool mmap) {
  if(binary) return readGraphFromBinary<vertex>(iFile,symmetric);
  else return readGraphFromFile<vertex>(iFile,symmetric,mmap);
}

template <class vertex>
graph<vertex> readCompressedGraph(char* fname, bool isSymmetric, bool mmap) {
  char* s;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    // Cannot mutate graph unless we copy.
    char *bytes = newA(char, S.n);
    parallel_for(size_t i=0; i<S.n; i++) {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    s = bytes;
  } else {
    ifstream in(fname,ifstream::in |ios::binary);
    in.seekg(0,ios::end);
    long size = in.tellg();
    in.seekg(0);
    cout << "size = " << size << endl;
    s = (char*) malloc(size);
    in.read(s,size);
    in.close();
  }

  long* sizes = (long*) s;
  long n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  cout << "n = "<<n<<" m = "<<m<<" totalSpace = "<<totalSpace<<endl;
  cout << "reading file..."<<endl;

  uintT* offsets = (uintT*) (s+3*sizeof(long));
  long skip = 3*sizeof(long) + (n+1)*sizeof(intT);
  uintE* Degrees = (uintE*) (s+skip);
  skip+= n*sizeof(intE);
  uchar* edges = (uchar*)(s+skip);

  uintT* inOffsets;
  uchar* inEdges;
  uintE* inDegrees;
  if(!isSymmetric){
    skip += totalSpace;
    uchar* inData = (uchar*)(s + skip);
    sizes = (long*) inData;
    long inTotalSpace = sizes[0];
    cout << "inTotalSpace = "<<inTotalSpace<<endl;
    skip += sizeof(long);
    inOffsets = (uintT*) (s + skip);
    skip += (n+1)*sizeof(uintT);
    inDegrees = (uintE*)(s+skip);
    skip += n*sizeof(uintE);
    inEdges = (uchar*)(s + skip);
  } else {
    inOffsets = offsets;
    inEdges = edges;
    inDegrees = Degrees;
  }


  vertex *V = newA(vertex,n);
  parallel_for(long i=0;i<n;i++) {
    long o = offsets[i];
    uintT d = Degrees[i];
    V[i].setOutDegree(d);
    V[i].setOutNeighbors(edges+o);
  }

  if(sizeof(vertex) == sizeof(compressedAsymmetricVertex)){
    parallel_for(long i=0;i<n;i++) {
      long o = inOffsets[i];
      uintT d = inDegrees[i];
      V[i].setInDegree(d);
      V[i].setInNeighbors(inEdges+o);
    }
  }

  cout << "creating graph..."<<endl;
  Compressed_Mem<vertex>* mem = new Compressed_Mem<vertex>(V, s);

  graph<vertex> G(V,n,m,mem);
  return G;
}

/* prefix sum used by the preprocess function defined below */
static pvector<uintT> ParallelPrefixSum (const pvector<uintT> &degrees) {
  const size_t block_size = 1<<20;
  const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
  pvector<uintT> local_sums(num_blocks);
  #pragma omp parallel for
  for (size_t block=0; block < num_blocks; block++) {
    uintT lsum = 0;
    size_t block_end = std::min((block + 1) * block_size, degrees.size());
    for (size_t i=block * block_size; i < block_end; i++)
      lsum += degrees[i];
    local_sums[block] = lsum;
  }
  pvector<uintT> bulk_prefix(num_blocks+1);
  uintT total = 0;
  for (size_t block=0; block < num_blocks; block++) {
    bulk_prefix[block] = total;
    total += local_sums[block];
  }
  bulk_prefix[num_blocks] = total;
  pvector<uintT> prefix(degrees.size() + 1);
  #pragma omp parallel for
  for (size_t block=0; block < num_blocks; block++) {
    uintT local_total = bulk_prefix[block];
    size_t block_end = std::min((block + 1) * block_size, degrees.size());
    for (size_t i=block * block_size; i < block_end; i++) {
      prefix[i] = local_total;
      local_total += degrees[i];
    }
  }
  prefix[degrees.size()] = bulk_prefix[num_blocks];
  return prefix;
}

template <class vertex>
void generateHubClusterMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto numVertices = GA.n;
    auto numEdges    = GA.m;
    vertex *origG    = GA.V;
    pvector<degree_nodeid_t> degree_id_pairs(numVertices);
    uintT avgDegree  = numEdges / numVertices;
    uintT hubCount {0};

    const int PADDING = 64 / sizeof(uintE);
    uintE* localOffsets = new uintE[omp_get_max_threads() * PADDING]();
    uintE partitionSz = numVertices / omp_get_max_threads();

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uintE startID = partitionSz * tid;
        uintE stopID  = partitionSz * (tid + 1);
        if (tid == omp_get_max_threads() - 1) {
            stopID = numVertices;    
        }
        for (uintE n = startID; n < stopID; ++n) {
            vertex vtx = origG[n];
            if (useOutdeg) {
                if (vtx.getOutDegree() > avgDegree) {
                    ++localOffsets[tid * PADDING];
                    new_ids[n] = 1;
                }
            }
            else {
                if (vtx.getInDegree() > avgDegree) {
                    ++localOffsets[tid * PADDING];
                    new_ids[n] = 1;
                }
            }
        }
    }
    uintE sum {0};
    for (int tid = 0; tid < omp_get_max_threads(); ++tid) {
        auto origCount = localOffsets[tid * PADDING];
        localOffsets[tid * PADDING] = sum;
        sum += origCount;
    }

    /* Step II - assign a remap for the hub vertices first */
#pragma omp parallel
    {
        uintE localCtr {0};
        int tid = omp_get_thread_num();
        uintE startID = partitionSz * tid;
        uintE stopID  = partitionSz * (tid + 1);
        if (tid == omp_get_max_threads() - 1) {
            stopID = numVertices;
        }
        for (uintE n = startID; n < stopID; ++n) {
            if (new_ids[n] != UINT_E_MAX) {
                new_ids[n] = localOffsets[tid * PADDING] + localCtr;
                ++localCtr;
            }
        }
    }
    delete[] localOffsets;

    /* Step III - assigning a remap for (easy) non hub vertices */
    auto numHubs = sum;
    SlidingQueue<uintE> queue(numHubs);
#pragma omp parallel
    {
        //assert(omp_get_max_threads() == 56);
        QueueBuffer<uintE> lqueue(queue, numHubs / omp_get_max_threads());
#pragma omp for
        for (uintE n = numHubs; n < numVertices; ++n) {
            if (new_ids[n] == UINT_E_MAX) {
                // This steps preserves the ordering of the original graph (as much as possible)
                new_ids[n] = n;
            }
            else {
                uintE remappedTo = new_ids[n];
                if (new_ids[remappedTo] == UINT_E_MAX) {
                    // safe to swap Ids because the original vertex is a non-hub 
                    new_ids[remappedTo] = n;
                }
                else {
                    // Cannot swap ids because original vertex was a hub (swapping 
                    // would disturb sorted ordering of hubs - not allowed)
                    lqueue.push_back(n);
                }
            }
        }
        lqueue.flush(); 
    }
    queue.slide_window(); //the queue keeps a list of vertices where a simple swap of locations is not possible

    /* Step IV - assigning remaps for remaining non hubs */
    uintE unassignedCtr {0};
    auto q_iter = queue.begin();
#pragma omp parallel for
    for (uintE n = 0; n < numHubs; ++n) {
        if (new_ids[n] == UINT_E_MAX) {
            uintE u = *(q_iter + __sync_fetch_and_add(&unassignedCtr, 1));
            new_ids[n] = u;
        }
    }

    t.Stop();
    t.PrintTime("HubCluster Map Time", t.Seconds());
}


template <class vertex>
void generateHubSortMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto numVertices = GA.n;
    auto numEdges    = GA.m;
    vertex *origG    = GA.V;
    pvector<degree_nodeid_t> degree_id_pairs(numVertices);
    uintT avgDegree  = numEdges / numVertices;
    uintT hubCount {0};

   /* STEP I - collect degrees of all vertices */
#pragma omp parallel for reduction(+ : hubCount)
    for (uintE v = 0; v < numVertices; ++v) {
        vertex vtx = origG[v];
        if (useOutdeg) {
            degree_id_pairs[v] = std::make_pair(vtx.getOutDegree(), v);
            if (vtx.getOutDegree() > avgDegree) {
                ++hubCount;
            }
        }
        else {
            degree_id_pairs[v] = std::make_pair(vtx.getInDegree(), v);
            if (vtx.getInDegree() > avgDegree) {
                ++hubCount;
            }
        }
    }

    /* Step II - sort the degrees in parallel */
    __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
            std::greater<degree_nodeid_t>());

    /* Step III - make a remap based on the sorted degree list [Only for hubs] */
#pragma omp parallel for
    for (uintE n = 0; n < hubCount; ++n) {
        new_ids[degree_id_pairs[n].second] = n;
    }
    //clearing space from degree pairs
    pvector<degree_nodeid_t>().swap(degree_id_pairs);

    /* Step IV - assigning a remap for (easy) non hub vertices */
    auto numHubs = hubCount;
    SlidingQueue<uintE> queue(numHubs);
#pragma omp parallel
    {
        QueueBuffer<uintE> lqueue(queue, numHubs / omp_get_max_threads());
#pragma omp for
        for (uintE n = numHubs; n < numVertices; ++n) {
            if (new_ids[n] == UINT_E_MAX) {
                // This steps preserves the ordering of the original graph (as much as possible)
                new_ids[n] = n;
            }
            else {
                uintE remappedTo = new_ids[n];
                if (new_ids[remappedTo] == UINT_E_MAX) {
                    // safe to swap Ids because the original vertex is a non-hub 
                    new_ids[remappedTo] = n;
                }
                else {
                    // Cannot swap ids because original vertex was a hub (swapping 
                    // would disturb sorted ordering of hubs - not allowed)
                    lqueue.push_back(n);
                }
            }
        }
        lqueue.flush(); 
    }
    queue.slide_window(); //the queue keeps a list of vertices where a simple swap of locations is not possible
    /* Step V - assigning remaps for remaining non hubs */
    uintE unassignedCtr {0};
    auto q_iter = queue.begin();
#pragma omp parallel for
    for (uintE n = 0; n < numHubs; ++n) {
        if (new_ids[n] == UINT_E_MAX) {
            uintE u = *(q_iter + __sync_fetch_and_add(&unassignedCtr, 1));
            new_ids[n] = u;
        }
    }

    t.Stop();
    t.PrintTime("HubSort Map Time", t.Seconds());
}


/* 
  Preprocess a graph based on outdegrees or indegrees
  PageRank-specific Optimizations for directed graphs - 
  1) We do not create a new outNeighbors list (because it pull-only)
  2) We only create new out-degrees because PR uses it during computation
*/
template <class vertex>
graph<vertex> preprocessGraph(graph<vertex> GA, bool isSym, bool useOutdeg, 
                              pvector<uintE>& new_ids, bool isPageRank = false, 
                              bool isDenseWrite = false, ReorderingAlgo reordering_algo = DBG)
{
    Timer t; 
    t.Start();
    auto numVertices = GA.n;
    auto numEdges    = GA.m;
    vertex *origG    = GA.V;
    
    if (!isSym) {

        generateMapping(GA, reordering_algo, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);

        /* Step VI - generate degree list for new graph */
        pvector<uintT> degrees(numVertices);
        pvector<uintT> inv_degrees(numVertices);
        #pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            auto newID = new_ids[v]; 
            if (useOutdeg) {
                vertex vtx         = origG[v];
                degrees[newID]     = vtx.getOutDegree(); 
                inv_degrees[newID] = vtx.getInDegree();
            }
            else {
                vertex vtx         = origG[v];
                degrees[newID]     = vtx.getInDegree(); 
                inv_degrees[newID] = vtx.getOutDegree();
            }
        }

        /* Step VII - make a new vertex list for the new graph */
		pvector<uintT> offsets     = ParallelPrefixSum(degrees);
		pvector<uintT> inv_offsets = ParallelPrefixSum(inv_degrees);
        //clearing space from degree lists
        pvector<uintT>().swap(degrees);
        pvector<uintT>().swap(inv_degrees);
        #ifndef WEIGHTED
            uintE* outEdges = newA(uintE, numEdges);
            uintE* inEdges  = newA(uintE, numEdges);
        #else
            intE* outEdges = newA(intE, 2 * numEdges);
            intE* inEdges  = newA(intE, 2 * numEdges);
        #endif
        vertex* newV    = newA(vertex, numVertices);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (uintE v = 0; v < numVertices; ++v) {
            /* note that vertex IDs u and v belong to the space of original vertex IDs */
            if (!isPageRank) {
                //copy out-neighbors
                auto newID = new_ids[v];
                newV[newID].setOutDegree(origG[v].getOutDegree());
                #ifndef WEIGHTED
                    if (useOutdeg)
                        newV[newID].setOutNeighbors(outEdges + offsets[newID]); 
                    else
                        newV[newID].setOutNeighbors(outEdges + inv_offsets[newID]); 
                #else
                    if (useOutdeg)
                        newV[newID].setOutNeighbors(outEdges + 2 * offsets[newID]); 
                    else
                        newV[newID].setOutNeighbors(outEdges + 2 * inv_offsets[newID]); 
                #endif
                for (uintE u = 0; u < origG[v].getOutDegree(); ++u) {
                    auto origNgh = origG[v].getOutNeighbor(u);
                    newV[newID].setOutNeighbor(u, new_ids[origNgh]);
                    #ifdef WEIGHTED
                        newV[newID].setOutWeight(u, origG[v].getOutWeight(u));
                    #endif
                }
                if (!isDenseWrite) {
                    /* for dense-write pushonly apps we dont need in-neighbors */
                    //copy in-neighbors
                    newV[newID].setInDegree(origG[v].getInDegree());
                    #ifndef WEIGHTED
                        if (useOutdeg) 
                            newV[newID].setInNeighbors(inEdges + inv_offsets[newID]); 
                        else
                            newV[newID].setInNeighbors(inEdges + offsets[newID]); 
                    #else
                        if (useOutdeg) 
                            newV[newID].setInNeighbors(inEdges + 2 * inv_offsets[newID]); 
                        else
                            newV[newID].setInNeighbors(inEdges + 2 * offsets[newID]); 
                    #endif
                    for (uintE u = 0; u < origG[v].getInDegree(); ++u) {
                        auto origNgh = origG[v].getInNeighbor(u);
                        newV[newID].setInNeighbor(u, new_ids[origNgh]);
                        #ifdef WEIGHTED
                            newV[newID].setInWeight(u, origG[v].getInWeight(u));
                        #endif
                    }
                }
            }
            else {
                /* PageRank - no need to apply weighted conditionals */
                //copy in-neighbors
                auto newID = new_ids[v];
                newV[newID].setInDegree(origG[v].getInDegree());
                if (useOutdeg) 
                    newV[newID].setInNeighbors(inEdges + inv_offsets[newID]); 
                else
                    newV[newID].setInNeighbors(inEdges + offsets[newID]); 
                for (uintE u = 0; u < origG[v].getInDegree(); ++u) {
                    auto origNgh = origG[v].getInNeighbor(u);
                    newV[newID].setInNeighbor(u, new_ids[origNgh]);
                }
                //only set out-degrees
                newV[newID].setOutDegree(origG[v].getOutDegree());
            }
        }

        /* Step V - make the new graph */ 
        Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(newV,numVertices,numEdges,outEdges,inEdges);
        t.Stop();
        string s = ReorderingAlgoStr(reordering_algo) + " Total Map Time";
        t.PrintTime(s.c_str(), t.Seconds());
        return graph<vertex>(newV,numVertices,numEdges,mem);
    }
    else {
        /* undirected graph */

        generateMapping(GA, reordering_algo, isSym, true, new_ids, isPageRank, isDenseWrite);

        /* Step VI - generate degree list for new graph */
        pvector<uintT> degrees(numVertices);
        #pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            auto newID = new_ids[v]; 
            vertex vtx         = origG[v];
            degrees[newID]     = vtx.getOutDegree(); 
        }

        /* Step VII - make a new vertex list for the new graph */
		pvector<uintT> offsets     = ParallelPrefixSum(degrees);
        //clearing space from degrees
        pvector<uintT>().swap(degrees);
        #ifndef WEIGHTED
            uintE* outEdges = newA(uintE, numEdges);
        #else
            intE* outEdges = newA(intE, 2 * numEdges);
        #endif
        vertex* newV    = newA(vertex, numVertices);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (uintE v = 0; v < numVertices; ++v) {
            /* note that vertex IDs u and v belong to the space of original vertex IDs */
            //copy neighbors
            auto newID = new_ids[v];
            newV[newID].setOutDegree(origG[v].getOutDegree());
            #ifndef WEIGHTED
                newV[newID].setOutNeighbors(outEdges + offsets[newID]); 
            #else
                newV[newID].setOutNeighbors(outEdges + 2 * offsets[newID]); 
            #endif
            for (uintE u = 0; u < origG[v].getOutDegree(); ++u) {
                auto origNgh = origG[v].getOutNeighbor(u);
                newV[newID].setOutNeighbor(u, new_ids[origNgh]);
                #ifdef WEIGHTED
                    newV[newID].setOutWeight(u, origG[v].getOutWeight(u));
                #endif
            }
        }

        /* Step V - make the new graph */ 
        Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(newV,numVertices,numEdges,outEdges);
        t.Stop();
        string s = ReorderingAlgoStr(reordering_algo) + " Total Map Time";
        t.PrintTime(s.c_str(), t.Seconds());
        return graph<vertex>(newV,numVertices,numEdges,mem);
    }
}

