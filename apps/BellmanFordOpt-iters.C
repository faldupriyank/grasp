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
#define WEIGHTED 1
#include "ligra.h"

// Added by Priyank to increase spatial locality
struct BF_Data {
    intE ShortestPathLen;
    int Visited;
};
struct BF_F {
  BF_Data* Packed_Len_Visited;
  intE* ShortestPathLen;
  int* Visited;
  BF_F(BF_Data* _Packed_Len_Visited) :
    Packed_Len_Visited(_Packed_Len_Visited) {}
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update ShortestPathLen if found a shorter path
    intE newDist = Packed_Len_Visited[s].ShortestPathLen + edgeLen;
    if(Packed_Len_Visited[d].ShortestPathLen > newDist) {
      Packed_Len_Visited[d].ShortestPathLen = newDist;
      if(Packed_Len_Visited[d].Visited == 0) { Packed_Len_Visited[d].Visited = 1 ; return 1;}
    }
    return 0;
  }
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic Update
    intE newDist = Packed_Len_Visited[s].ShortestPathLen + edgeLen;
    return (writeMin(&Packed_Len_Visited[d].ShortestPathLen,newDist) &&
	    CAS(&Packed_Len_Visited[d].Visited,0,1));
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

//reset visited vertices
struct BF_Vertex_F {
  BF_Data* Packed_Len_Visited;
  BF_Vertex_F(BF_Data* _Packed_Len_Visited) : Packed_Len_Visited(_Packed_Len_Visited) {}
  inline bool operator() (uintE i){
    Packed_Len_Visited[i].Visited = 0;
    return 1;
  }
};

#ifdef _OUTPUT_
void sampleOutput(BF_Data* arr, bool preprocessed, const pvector<uintE> &new_ids) {
    std::cout << "BellmanFord output! data preprocessed? " << preprocessed << std::endl;
    for ( long i = 0 ; i < 10 && i < new_ids.size() ; i++ ) {
        if ( preprocessed )
            std::cout << i << " " << arr[new_ids[i]].ShortestPathLen << std::endl;
        else
            std::cout << i << " " << arr[i].ShortestPathLen << std::endl;
    }
}
#endif

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P, pvector<uintE> &new_ids) {
#ifdef _SNIPER_
    std::cout << "_SNIPER_" << std::endl;
#endif
  Timer tm;
  tm.Start();
  ulong sourceCtr {0};
  const long NUM_ITERS = P.getOptionLongValue("-num_roots", 8);
  for (int iteration = 0; iteration < NUM_ITERS; ++iteration) {
      bool preprocessed = (new_ids[0] != new_ids[1]); 
      long n = GA.n;
      long start;
      if (preprocessed) {
          do {
              start = hashInt(sourceCtr++) % n; //random numbers are in original vertex space
          } while (GA.V[new_ids[start]].getOutDegree() == 0);
      }
      else {
          do {
              start = hashInt(sourceCtr++) % n;
          } while (GA.V[start].getOutDegree() == 0);
      }

      if (preprocessed)
          start = new_ids[start];
  //initialize ShortestPathLen to "infinity"
#ifndef ALIGNED
      BF_Data* Packed_Len_Visited = newA(BF_Data,n);
#else
      BF_Data* Packed_Len_Visited {nullptr};
      posix_memalign((void**)&Packed_Len_Visited, 64, n * sizeof(BF_Data));
      assert(Packed_Len_Visited != nullptr && ((uintptr_t)Packed_Len_Visited % 64 == 0) && "App Malloc Failure\n");
#endif
      {parallel_for(long i=0;i<n;i++) {
                                          Packed_Len_Visited[i].ShortestPathLen = INT_MAX/2;
                                          Packed_Len_Visited[i].Visited = 0;
                                      }}
      Packed_Len_Visited[start].ShortestPathLen = 0;


      vertexSubset Frontier(n,start); //initial frontier

#ifdef _SNIPER_
  // No need for frontier_frac as it is push-based only model
  // as frontier is not acessed using irregular accesses.
  if ( dense_iter_nooutput == 0 && dense_iter_nooutput == 0 ) {
      single_roi_output = 1;
      single_roi_nooutput = 1;
      add_region("propertyA", Packed_Len_Visited, 100, n);
  }
#endif


      long round = 0;
      while(!Frontier.isEmpty()){
          if(round == n) {
              //negative weight cycle
              {parallel_for(long i=0;i<n;i++) Packed_Len_Visited[i].ShortestPathLen = -(INT_E_MAX/2);}
              break;
          }
          vertexSubset output = edgeMap(GA, Frontier, BF_F(Packed_Len_Visited), GA.m/20, dense_forward);
          vertexMap(output,BF_Vertex_F(Packed_Len_Visited));
          Frontier.del();
          Frontier = output;
          round++;
      }
      Frontier.del();
      free(Packed_Len_Visited);
  }
  tm.Stop();
  tm.PrintTime("BellmanFord-iters Run Time(sec) ", tm.Seconds()); 
}
