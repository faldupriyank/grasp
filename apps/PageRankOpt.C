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
#include "ligra.h"
#include "math.h"

// Added by Priyank to increase spatial locality
struct PR_Data {
    double p_curr;
    uintT degree;
};

template <class vertex>
struct PR_F {
  double* p_next;
  PR_Data* Packed_Curr_Degree;
  PR_F(PR_Data* _Packed_Curr_Degree, double* _p_next) : 
    Packed_Curr_Degree(_Packed_Curr_Degree), p_next(_p_next) {}
  inline bool update(uintE s, uintE d){ //update function applies PageRank equation
    p_next[d] += Packed_Curr_Degree[s].p_curr/Packed_Curr_Degree[s].degree;
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    writeAdd(&p_next[d],Packed_Curr_Degree[s].p_curr/Packed_Curr_Degree[s].degree);
    return 1;
  }
  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
  double damping;
  double addedConstant;
  double* p_next;
  PR_Vertex_F(double* _p_next, double _damping, intE n) :
    p_next(_p_next), 
    damping(_damping), addedConstant((1-_damping)*(1/(double)n)){}
  inline bool operator () (uintE i) {
    p_next[i] = damping*p_next[i] + addedConstant;
    return 1;
  }
};

struct PR_Vertex_Swap {
  PR_Data* Packed_Curr_Degree;
  double* p_next;
  PR_Vertex_Swap(PR_Data* _Packed_Curr_Degree, double* _p_next) :
    Packed_Curr_Degree(_Packed_Curr_Degree), p_next(_p_next) {}
  inline bool operator () (uintE i) {
    Packed_Curr_Degree[i].p_curr = p_next[i];
    p_next[i] = 0.0;
    return 1;
  }
};

#ifdef _OUTPUT_

void sampleOutput(double* arr, long numElements, pvector<uintE> &new_ids) {
    bool preprocessed = (new_ids[0] != new_ids[1]);
    std::cout << "PageRank output! data preprocessed? " << preprocessed << std::endl;
    if ( !preprocessed ) {
        for (long v = 0; v < 10 && v < numElements; ++v) {
            std::cout << v << " " << arr[v] << "\n";
        }
    }
    else {
        for (long v = 0; v < 10 && v < numElements; ++v) {
            std::cout << v << " " << arr[new_ids[v]] << "\n";
        }
    }
    return;

}
void sampleOutput(PR_Data* arr, long numElements, pvector<uintE> &new_ids) {
    bool preprocessed = (new_ids[0] != new_ids[1]);
    std::cout << "PageRank output! data preprocessed? " << preprocessed << std::endl;
    if ( !preprocessed ) {
        for (long v = 0; v < 10 && v < numElements; ++v) {
            std::cout << v << " " << arr[v].p_curr << "\n";
        }
    }
    else {
        for (long v = 0; v < 10 && v < numElements; ++v) {
            std::cout << v << " " << arr[new_ids[v]].p_curr << "\n";
        }
    }
    return;
}
#endif

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P, pvector<uintE> &new_ids) {
#ifdef _SNIPER_
    std::cout << "_SNIPER_" << std::endl;
#endif
  Timer tm;
  tm.Start();
  long maxIters = P.getOptionLongValue("-maxiters",100);
  const intE n = GA.n;
  const double damping = 0.85, epsilon = 0.0000001;
  
  double one_over_n = 1/(double)n;
  #ifndef ALIGNED
  PR_Data* Packed_Curr_Degree = newA(PR_Data,n);
  #else
  double* Packed_Curr_Degree {nullptr};
  posix_memalign((void**) &Packed_Curr_Degree, 64, sizeof(double) * n);
  assert(Packed_Curr_Degree != nullptr && ((uintptr_t)Packed_Curr_Degree % 64 == 0) && "App Malloc Failure\n");
  #endif
  double* p_next = newA(double,n);
  {parallel_for(long i=0;i<n;i++) {
      Packed_Curr_Degree[i].p_curr = one_over_n;
      Packed_Curr_Degree[i].degree = GA.V[i].getOutDegree();
      p_next[i] = 0;
  }} //0 if unchanged
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}

  vertexSubset Frontier(n,n,frontier);
  
#ifdef _SNIPER_
  single_roi_output = 1;
  single_roi_nooutput = 1;
  add_region("propertyA", Packed_Curr_Degree, 100-frontier_frac, n);
#endif
  long iter = 0;
  while(iter++ < maxIters) {
    edgeMap(GA,Frontier,PR_F<vertex>(Packed_Curr_Degree,p_next),0, no_output);
    vertexMap(Frontier,PR_Vertex_F(p_next,damping,n));
    //compute L1-norm between p_curr and p_next
    double L1_norm = 0;

#pragma omp parallel for reduction(+:L1_norm)
    for(long i=0;i<n;i++) {
      Packed_Curr_Degree[i].p_curr = fabs(Packed_Curr_Degree[i].p_curr-p_next[i]);
      L1_norm += Packed_Curr_Degree[i].p_curr;
    }
    if(L1_norm < epsilon) break;
    //reset p_curr
    vertexMap(Frontier,PR_Vertex_Swap(Packed_Curr_Degree, p_next));
  }
  Frontier.del();
#ifdef _OUTPUT_
  tm.Stop();
  if ( iter < maxIters ) {
      sampleOutput(p_next, n, new_ids);
  } else {
      sampleOutput(Packed_Curr_Degree, n, new_ids);
  }
  free(Packed_Curr_Degree); free(p_next); 
#else
  free(Packed_Curr_Degree); free(p_next);
  tm.Stop();
#endif
  tm.PrintTime("PageRank Run Time(sec) ", tm.Seconds()); 
}
