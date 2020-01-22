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
    double delta;
    uintT degree;
};
template <class vertex>
struct PR_F {
  PR_Data* Packed_Delta_Degree;
  double *nghSum;

  PR_F(PR_Data* _Packed_Delta_Degree, double* _nghSum) : 
      Packed_Delta_Degree(_Packed_Delta_Degree), nghSum(_nghSum) {}

  inline bool update(uintE s, uintE d){
    double oldVal = nghSum[d];
    nghSum[d] += (Packed_Delta_Degree[s].delta/Packed_Delta_Degree[s].degree);
    return oldVal == 0;
  }
  inline bool updateAtomic (uintE s, uintE d) {
    volatile double oldV, newV; 
    do { //basically a fetch-and-add
      oldV = nghSum[d]; newV = oldV + (Packed_Delta_Degree[s].delta/Packed_Delta_Degree[s].degree);
    } while(!CAS(&nghSum[d],oldV,newV));
    return oldV == 0.0;
  }
  inline bool cond (uintE d) { return cond_true(d); }};

struct PR_Vertex_F_FirstRound {
  double damping, addedConstant, one_over_n, epsilon2;
  double* p, *nghSum;
  PR_Data* Packed_Delta_Degree;
  PR_Vertex_F_FirstRound(double* _p, PR_Data* _Packed_Delta_Degree, double* _nghSum, double _damping, double _one_over_n,double _epsilon2) :
    p(_p),
    damping(_damping), Packed_Delta_Degree(_Packed_Delta_Degree), nghSum(_nghSum), one_over_n(_one_over_n),
    addedConstant((1-_damping)*_one_over_n),
    epsilon2(_epsilon2) {}
  inline bool operator () (uintE i) {
    Packed_Delta_Degree[i].delta = damping*nghSum[i]+addedConstant;
    p[i] += Packed_Delta_Degree[i].delta;
    Packed_Delta_Degree[i].delta -=one_over_n; //subtract off delta from initialization
    return (fabs(Packed_Delta_Degree[i].delta) > epsilon2 * p[i]);
  }
};
void sampleOutput(double* arr, long numElements, pvector<uintE> &new_ids) {
    bool preprocessed = (new_ids[0] != new_ids[1]);
    std::cout << "PageRankDelta output! data preprocessed? " << preprocessed << std::endl;
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
struct PR_Vertex_F {
  double damping, epsilon2;
  double* p, *nghSum;
  PR_Data* Packed_Delta_Degree;
  PR_Vertex_F(double* _p, PR_Data* _Packed_Delta_Degree, double* _nghSum, double _damping, double _epsilon2) :
    p(_p),
    damping(_damping), Packed_Delta_Degree(_Packed_Delta_Degree), nghSum(_nghSum), 
    epsilon2(_epsilon2) {}
  inline bool operator () (uintE i) {
    Packed_Delta_Degree[i].delta = nghSum[i]*damping;
    if (fabs(Packed_Delta_Degree[i].delta) > epsilon2*p[i]) { p[i]+=Packed_Delta_Degree[i].delta; return 1;}
    else return 0;
  }
};

struct PR_Vertex_Reset {
  double* nghSum;
  PR_Vertex_Reset(double* _nghSum) :
    nghSum(_nghSum) {}
  inline bool operator () (uintE i) {
    nghSum[i] = 0.0;
    return 1;
  }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P, pvector<uintE> &new_ids) {
#ifdef _SNIPER_
    std::cout << "_SNIPER_" << std::endl;
#endif
  Timer tm;
  tm.Start();
  long maxIters = P.getOptionLongValue("-maxiters",100);
  const long n = GA.n;
  const double damping = 0.85;
  const double epsilon = 0.0000001;
  const double epsilon2 = 0.01;

  double one_over_n = 1/(double)n;
  #ifndef ALIGNED
  double* p = newA(double,n), *nghSum = newA(double,n);
  PR_Data *Packed_Delta_Degree= newA(PR_Data,n);
  #else
  double* p = newA(double,n), *nghSum = newA(double,n);
  PR_Data* Packed_Delta_Degree {null};
  posix_memalign((void**) &Packed_Delta_Degree, 64, sizeof(PR_Data) * n);
  assert(Packed_Delta_Degree!= nullptr && ((uintptr_t)Packed_Delta_Degree % 64 == 0) && "App Malloc Failure\n");
  #endif
  bool* frontier = newA(bool,n);
  parallel_for(long i=0;i<n;i++) {
    p[i] = 0.0;//one_over_n;
    Packed_Delta_Degree[i].delta = one_over_n; //initial delta propagation from each vertex
    Packed_Delta_Degree[i].degree = GA.V[i].getOutDegree();
    nghSum[i] = 0.0;
    frontier[i] = 1;
  }

  vertexSubset Frontier(n,n,frontier);
  bool* all = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) all[i] = 1;}
  vertexSubset All(n,n,all); //all vertices

#ifdef _SNIPER_
  single_roi_output = 1;
  single_roi_nooutput = 1;
  add_region("propertyA", Packed_Delta_Degree, 100-frontier_frac, n);
#endif
  long round = 0;
  while(round++ < maxIters) {
    edgeMap(GA,Frontier,PR_F<vertex>(Packed_Delta_Degree,nghSum),GA.m/20, no_output);
    vertexSubset active 
      = (round == 1) ? 
      vertexFilter(All,PR_Vertex_F_FirstRound(p,Packed_Delta_Degree,nghSum,damping,one_over_n,epsilon2)) :
      vertexFilter(All,PR_Vertex_F(p,Packed_Delta_Degree,nghSum,damping,epsilon2));
    //compute L1-norm (use nghSum as temp array)
    {parallel_for(long i=0;i<n;i++) {
      nghSum[i] = fabs(Packed_Delta_Degree[i].delta); }}
    double L1_norm = sequence::plusReduce(nghSum,n);
    if(L1_norm < epsilon) break;
    //reset
    vertexMap(All,PR_Vertex_Reset(nghSum));
    Frontier.del();
    Frontier = active;
  }
  Frontier.del(); free(Packed_Delta_Degree); free(nghSum); All.del();
#ifdef _OUTPUT_
  tm.Stop();
  sampleOutput(p, n, new_ids);
  free(p); 
#else
  free(p); 
  tm.Stop();
#endif
  tm.PrintTime("PageRankDelta Run Time(sec) ", tm.Seconds());
}
