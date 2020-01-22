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

struct BF_F {
  intE* ShortestPathLen;
  int* Visited;
  BF_F(intE* _ShortestPathLen, int* _Visited) : 
    ShortestPathLen(_ShortestPathLen), Visited(_Visited) {}
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update ShortestPathLen if found a shorter path
    intE newDist = ShortestPathLen[s] + edgeLen;
    if(ShortestPathLen[d] > newDist) {
      ShortestPathLen[d] = newDist;
      if(Visited[d] == 0) { Visited[d] = 1 ; return 1;}
    }
    return 0;
  }
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic Update
    intE newDist = ShortestPathLen[s] + edgeLen;
    return (writeMin(&ShortestPathLen[d],newDist) &&
	    CAS(&Visited[d],0,1));
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

//reset visited vertices
struct BF_Vertex_F {
  int* Visited;
  BF_Vertex_F(int* _Visited) : Visited(_Visited) {}
  inline bool operator() (uintE i){
    Visited[i] = 0;
    return 1;
  }
};


template <class vertex>
void Compute(graph<vertex>& GA, commandLine P, pvector<uintE> &new_ids) {
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
      intE* ShortestPathLen = newA(intE,n);
#else
      intE* ShortestPathLen {nullptr};
      posix_memalign((void**)&ShortestPathLen, 64, n * sizeof(intE));
      assert(ShortestPathLen != nullptr && ((uintptr_t)ShortestPathLen % 64 == 0) && "App Malloc Failure\n");
#endif
      {parallel_for(long i=0;i<n;i++) ShortestPathLen[i] = INT_MAX/2;}
      ShortestPathLen[start] = 0;

      int* Visited = newA(int,n);
      {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}

      vertexSubset Frontier(n,start); //initial frontier

      long round = 0;
      while(!Frontier.isEmpty()){
          if(round == n) {
              //negative weight cycle
              {parallel_for(long i=0;i<n;i++) ShortestPathLen[i] = -(INT_E_MAX/2);}
              break;
          }
          vertexSubset output = edgeMap(GA, Frontier, BF_F(ShortestPathLen,Visited), GA.m/20, dense_forward);
          vertexMap(output,BF_Vertex_F(Visited));
          Frontier.del();
          Frontier = output;
          round++;
      }
      Frontier.del(); free(Visited); free(ShortestPathLen);
  }
  tm.Stop();
  tm.PrintTime("BellmanFord-iters Run Time(sec) ", tm.Seconds()); 
}
