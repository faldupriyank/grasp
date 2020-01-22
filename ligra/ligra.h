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
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <cassert>
#include "parallel.h"
#include "gettime.h"
#include "timer.h"  //timer from GAP
#include "utils.h"
#include "vertex.h"
#include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"
#include "index_map.h"
#include "edgeMap_utils.h"
#ifdef _SNIPER_
#include "sim_api.h"
#endif
using namespace std;

#if 1 // Added by Priyank
int rand_gran = 1;
int num_roots = 8;

#ifdef _SNIPER_
// First iteration that has minimum %sniper_min_threshold% percentage of active vertices
// will be selected for simulation.
// For some applications, fraction of active vertices are very low throught the application,
// and thus sniper_min_threshold should be set very low.
long sniper_min_threshold;
long sniper_max_threshold;
int single_roi_output = 0;
int single_roi_nooutput = 0;
int dense_iter_output = 0;
int dense_iter_nooutput = 0;
int frontier_frac = 50; // May need to adjust this for different applications.
#endif

string map_file = "";

#ifdef _SNIPER_
void SimNamedMarkerCommon(uint64_t val, const string& str) {
    cout << str.c_str() << " " << val << endl;
    SimNamedMarker(val, str.c_str());
}
template <typename T>
inline void add_region(const string& region_name, T* base_ptr, uint64_t frac, long n) {
    const string region_0 = region_name + "-0";
    const string region_1 = region_name + "-1";
    const string region_n = region_name + "-n";
    const string region_f = region_name + "-f";
    SimNamedMarkerCommon(reinterpret_cast<uintptr_t>(base_ptr+0), region_0);
    SimNamedMarkerCommon(reinterpret_cast<uintptr_t>(base_ptr+1), region_1);
    SimNamedMarkerCommon(reinterpret_cast<uintptr_t>(base_ptr+n), region_n);
    SimNamedMarkerCommon(frac, region_f);
    std::cout << region_name << " " << std::hex << base_ptr+0 << " " << base_ptr+1 << " " << base_ptr+n << std::dec << " " << (uint64_t)(base_ptr+1) - (uint64_t)base_ptr << " " << (uint64_t)(base_ptr+n) - (uint64_t)base_ptr  << " " << ((uint64_t)(base_ptr+n) - (uint64_t)base_ptr) / ((uint64_t)(base_ptr+1) - (uint64_t)base_ptr) << " " << frac << " +++++ " << std::endl;
}
#endif
#endif

//*****START FRAMEWORK*****

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_no_filter = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
inline bool should_output(const flags& fl) { return !(fl & no_output); }
const int dynChunkSz = 64; //chunk size for openmp's dynamic scheduling

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDense(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA.n;
  vertex *G = GA.V;
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_gen<data>(next);

#ifdef _SNIPER_
    const long all_vert = vertexSubset.numVertices();
    const long active_vert = vertexSubset.size();
    const long perct_vert = active_vert * 100.0 / all_vert;
    std::cout << dense_iter_output++ << " out_iter: " << active_vert << " " << all_vert << " " << perct_vert << std::endl;
    if ( single_roi_output != 0 && perct_vert >= sniper_min_threshold ) {
        assert(vertexSubset.isDense);
        add_region<bool>("propertyB", (bool*) vertexSubset.frontierAddr(0), frontier_frac, n);
        SimNamedMarker(0, "grasp-enable");
        single_roi_output = 2;
        std::cout << "Enabling ROI: " << perct_vert << std::endl;
        SimRoiStart();
    }
#endif
    #pragma omp parallel for schedule (dynamic, dynChunkSz)
    for (long v=0; v<n; v++) {
      std::get<0>(next[v]) = 0;
      if (f.cond(v)) {
        G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
#ifdef _SNIPER_
    if ( single_roi_output == 2 ) {
        SimRoiEnd();
        SimNamedMarker(0, "grasp-disable");
        single_roi_output = 0;
        std::cout << "Disabling ROI." << std::endl;
    }
#endif
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_nooutput_gen<data>();
#ifdef _SNIPER_
    const long all_vert = vertexSubset.numVertices();
    const long active_vert = vertexSubset.size();
    const long perct_vert = active_vert * 100.0 / all_vert;
    std::cout << dense_iter_nooutput++ << " no_iter: " << active_vert << " " << all_vert << " " << perct_vert << std::endl;
    if ( single_roi_nooutput != 0 && perct_vert >= sniper_min_threshold && perct_vert <= sniper_max_threshold ) {
        assert(vertexSubset.isDense);
        add_region<bool>("propertyB", (bool*) vertexSubset.frontierAddr(0), frontier_frac, n);
        SimNamedMarker(0, "grasp-enable");
        single_roi_nooutput = 2;
        std::cout << "Enabling ROI: " << perct_vert << std::endl;
        SimRoiStart();
    }
#endif
    #pragma omp parallel for schedule (dynamic, dynChunkSz)
    for (long v=0; v<n; v++) {
      if (f.cond(v)) {
        G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
#ifdef _SNIPER_
    if ( single_roi_nooutput == 2 ) {
        SimRoiEnd();
        SimNamedMarker(0, "grasp-disable");
        single_roi_nooutput = 0;
        std::cout << "Disabling ROI." << std::endl;
    }
#endif
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDenseForward(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA.n;
  vertex *G = GA.V;
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_forward_gen<data>(next);
    parallel_for(long i=0;i<n;i++) { std::get<0>(next[i]) = 0; }
#ifdef _SNIPER_
    const long all_vert = vertexSubset.numVertices();
    const long active_vert = vertexSubset.size();
    const long perct_vert = active_vert * 100.0 / all_vert;
    std::cout << dense_iter_output++ << " fwd_iter: " << active_vert << " " << all_vert << " " << perct_vert << std::endl;
    if ( single_roi_output != 0 && perct_vert >= sniper_min_threshold && perct_vert <= sniper_max_threshold ) {
        assert(vertexSubset.isDense);
        add_region<bool>("propertyB", (bool*) vertexSubset.frontierAddr(0), 0, n);
        SimNamedMarker(0, "grasp-enable");
        assert(single_roi_output == 1);
        single_roi_output = 2;
        std::cout << "Enabling ROI: " << perct_vert << std::endl;
        SimRoiStart();
    }
#endif
    #pragma omp parallel for schedule (dynamic, dynChunkSz)
    for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        G[i].decodeOutNgh(i, f, g);
      }
    }
#ifdef _SNIPER_
    if ( single_roi_output == 2 ) {
        SimRoiEnd();
        SimNamedMarker(0, "grasp-disable");
        single_roi_output = 0;
        std::cout << "Disabling ROI." << std::endl;
    }
#endif

    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_forward_nooutput_gen<data>();
#ifdef _SNIPER_
    const long all_vert = vertexSubset.numVertices();
    const long active_vert = vertexSubset.size();
    const long perct_vert = active_vert * 100.0 / all_vert;
    std::cout << dense_iter_output++ << " fwd_iter: " << active_vert << " " << all_vert << " " << perct_vert << std::endl;
    if ( single_roi_output != 0 && perct_vert >= sniper_min_threshold && perct_vert <= sniper_max_threshold ) {
        assert(vertexSubset.isDense);
        add_region<bool>("propertyB", (bool*) vertexSubset.frontierAddr(0), 0, n);
        SimNamedMarker(0, "grasp-enable");
        assert(single_roi_output == 1);
        single_roi_output = 2;
        std::cout << "Enabling ROI: " << perct_vert << std::endl;
        SimRoiStart();
    }
#endif
    #pragma omp parallel for schedule (dynamic, dynChunkSz)
    for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        G[i].decodeOutNgh(i, f, g);
      }
    }
#ifdef _SNIPER_
    if ( single_roi_output == 2 ) {
        SimRoiEnd();
        SimNamedMarker(0, "graphinity-disable");
        single_roi_output = 0;
        std::cout << "Disabling ROI." << std::endl;
    }
#endif

    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse(graph<vertex>& GA, vertex* frontierVertices, VS& indices,
        uintT* degrees, uintT m, F &f, const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  S* outEdges;
  long outEdgeCount = 0;

  if (should_output(fl)) {
    uintT* offsets = degrees;
    outEdgeCount = sequence::plusScan(offsets, offsets, m);
    outEdges = newA(S, outEdgeCount);
    auto g = get_emsparse_gen<data>(outEdges);
    #pragma omp parallel for schedule (dynamic, dynChunkSz)
    for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i), o = offsets[i];
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, o, f, g);
    }
  } else {
    auto g = get_emsparse_nooutput_gen<data>();
    #pragma omp parallel for schedule (dynamic, dynChunkSz)
    for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i);
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, 0, f, g);
    }
  }

  if (should_output(fl)) {
    S* nextIndices = newA(S, outEdgeCount);
    if (fl & remove_duplicates) {
      if (GA.flags == NULL) {
        GA.flags = newA(uintE, n);
        parallel_for(long i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
      }
      auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(outEdges[i]); };
      remDuplicates(get_key, GA.flags, outEdgeCount, n);
    }
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
    free(outEdges);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  } else {
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter(graph<vertex>& GA,
    vertex* frontierVertices, VS& indices, uintT* offsets, uintT m, F& f,
    const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  long outEdgeCount = sequence::plusScan(offsets, offsets, m);
  S* outEdges = newA(S, outEdgeCount);

  auto g = get_emsparse_no_filter_gen<data>(outEdges);

  // binary-search into scan to map workers->chunks
  size_t b_size = 10000;
  size_t n_blocks = nblocks(outEdgeCount, b_size);

  uintE* cts = newA(uintE, n_blocks+1);
  size_t* block_offs = newA(size_t, n_blocks+1);

  auto offsets_m = make_in_imap<uintT>(m, [&] (size_t i) { return offsets[i]; });
  auto lt = [] (const uintT& l, const uintT& r) { return l < r; };
  parallel_for(size_t i=0; i<n_blocks; i++) {
    size_t s_val = i*b_size;
    block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
  }
  block_offs[n_blocks] = m;
  #pragma omp parallel for schedule (dynamic, dynChunkSz / 8)
  for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      // start and end are offsets in [m]
      size_t start = block_offs[i];
      size_t end = block_offs[i+1];
      uintT start_o = offsets[start];
      uintT k = start_o;
      for (size_t j=start; j<end; j++) {
        uintE v = indices.vtx(j);
        size_t num_in = frontierVertices[j].decodeOutNghSparseSeq(v, k, f, g);
        k += num_in;
      }
      cts[i] = (k - start_o);
    } else {
      cts[i] = 0;
    }
  }

  long outSize = sequence::plusScan(cts, cts, n_blocks);
  cts[n_blocks] = outSize;

  S* out = newA(S, outSize);

  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      size_t start = block_offs[i];
      size_t start_o = offsets[start];
      size_t out_off = cts[i];
      size_t block_size = cts[i+1] - out_off;
      for (size_t j=0; j<block_size; j++) {
        out[out_off + j] = outEdges[start_o + j];
      }
    }
  }
  free(outEdges); free(cts); free(block_offs);

  if (fl & remove_duplicates) {
    if (GA.flags == NULL) {
      GA.flags = newA(uintE, n);
      parallel_for(size_t i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
    }
    auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(out[i]); };
    remDuplicates(get_key, GA.flags, outSize, n);
    S* nextIndices = newA(S, outSize);
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
    free(out);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  }
  return vertexSubsetData<data>(n, outSize, out);
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapData(graph<vertex>& GA, VS &vs, F f,
    intT threshold = -1, const flags& fl=0) {
  long numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  if (numVertices != vs.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  if (vs.size() == 0) return vertexSubsetData<data>(numVertices);
  vs.toSparse();
  uintT* degrees = newA(uintT, m);
  vertex* frontierVertices = newA(vertex,m);
  {parallel_for (size_t i=0; i < m; i++) {
    uintE v_id = vs.vtx(i);
    vertex v = G[v_id];
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
  }}

  uintT outDegrees = sequence::plusReduce(degrees, m);
  if (outDegrees == 0) return vertexSubsetData<data>(numVertices);
  if (m + outDegrees > threshold) {
    vs.toDense();
    free(degrees); free(frontierVertices);
    return (fl & dense_forward) ?
      edgeMapDenseForward<data, vertex, VS, F>(GA, vs, f, fl) :
      edgeMapDense<data, vertex, VS, F>(GA, vs, f, fl);
  } else {
    auto vs_out =
      (should_output(fl) && fl & sparse_no_filter) ? // only call snof when we output
      edgeMapSparse_no_filter<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl) :
      edgeMapSparse<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl);
    free(degrees); free(frontierVertices);
    return vs_out;
  }
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class vertex, class VS, class F>
vertexSubset edgeMap(graph<vertex> GA, VS& vs, F f,
    intT threshold = -1, const flags& fl=0) {
  return edgeMapData<pbbs::empty>(GA, vs, f, threshold, fl);
}

/* General function to print stats about frontier size */
template <class VS>
void frontierStats(VS& vs, long numVertices, bool KCore = false) {
    if (KCore) {
        double percent = (static_cast<double>(vs.size()) / static_cast<double>(numVertices)) * 100;
        if (vs.dense()) {
            std::cout << "PULL iteration. Frontier size = " << percent << std::endl;
        }
        else { 
            std::cout << "PUSH iteration. Frontier size = " << percent << std::endl;
        }
    }
    return;
}

// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
// Weighted graphs are not yet supported, but this should be easy to do.
template <class vertex, class P>
vertexSubsetData<uintE> packEdges(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  using S = tuple<uintE, uintE>;
  vs.toSparse();
  vertex* G = GA.V; long m = vs.numNonzeros(); long n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto degrees = array_imap<uintT>(m);
  granular_for(i, 0, m, (m > 2000), {
    uintE v = vs.vtx(i);
    degrees[i] = G[v].getOutDegree();
  });
  long outEdgeCount = pbbs::scan_add(degrees, degrees);
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }

  bool* bits = newA(bool, outEdgeCount);
  uintE* tmp1 = newA(uintE, outEdgeCount);
  uintE* tmp2 = newA(uintE, outEdgeCount);
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
    }
  }
  free(bits); free(tmp1); free(tmp2);
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}

template <class vertex, class P>
vertexSubsetData<uintE> edgeMapFilter(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  vs.toSparse();
  if (fl & pack_edges) {
    return packEdges<vertex, P>(GA, vs, p, fl);
  }
  vertex* G = GA.V; long m = vs.numNonzeros(); long n = vs.numRows();
  using S = tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = G[v].countOutNgh(v, p);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = G[v].countOutNgh(v, p);
    }
  }
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}



//*****VERTEX FUNCTIONS*****

template <class F, class VS, typename std::enable_if<
  !std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i), V.vtxData(i));
    }
  }
}

template <class VS, class F, typename std::enable_if<
  std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i);
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i));
    }
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

template <class F>
vertexSubset vertexFilter2(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  {parallel_for(size_t i=0; i<m; i++) {
    uintE v = V.vtx(i);
    bits[i] = filter(v);
  }}
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}

template <class data, class F>
vertexSubset vertexFilter2(vertexSubsetData<data> V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  parallel_for(size_t i=0; i<m; i++) {
    auto t = V.vtxAndData(i);
    bits[i] = filter(std::get<0>(t), std::get<1>(t));
  }
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}

#if 1 // Added by Priyank
void mergeTwoPreprocessingIndices(pvector<uintE>& first, pvector<uintE>& second, pvector<uintE>& new_ids) {
    assert(first.size() == second.size());
    assert(first.size() == new_ids.size());

    long int numVertices = first.size();
    {parallel_for(long i = 0 ; i < numVertices ; i++ ) {
        new_ids[i] = second[first[i]];
    }}
}
#endif

//cond function that always returns true
inline bool cond_true (intT d) { return 1; }

template<class vertex>
void Compute(graph<vertex>&, commandLine, pvector<uintE> &new_ids);

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  bool isPageRank = (P.getOptionIntValue("-is_pagerank", -1) == 1);
  bool isDenseWrite = (P.getOptionIntValue("-is_dense_write", -1) == 1);
  long rounds = P.getOptionLongValue("-rounds",3);
#if 1 // Added by Priyank
#ifdef _SNIPER_
  rounds = 0;
  sniper_min_threshold = P.getOptionIntValue("-sniper-min", 40);
  sniper_max_threshold = P.getOptionIntValue("-sniper-max", 100);
#endif
  int degree_used_for_reordering = P.getOptionIntValue("-degree_used_for_reordering", -1);
  long threads = P.getOptionIntValue("-threads", omp_get_max_threads());
  if ( threads != omp_get_max_threads() ) {
      omp_set_num_threads(threads);
  }

  rand_gran = P.getOptionIntValue("-rand_gran", 1);
  num_roots = P.getOptionIntValue("-num_roots", 8);
  map_file = P.getOptionValue("-map_file", "");

  string degree_used_for_reordering_str = "none";
  enum ReorderingAlgo reordering_algo = (ReorderingAlgo) P.getOptionLongValue("-reordering_algo", DBG); 
  if ( degree_used_for_reordering == 0 || degree_used_for_reordering == 1 ) {
      if ( reordering_algo != ORIGINAL ) {
          degree_used_for_reordering_str = (degree_used_for_reordering == 0) ? "out-degree" : "in-degree";
      } else {
          degree_used_for_reordering_str = "none";
      }
  } else {
      reordering_algo = ORIGINAL;
  }
  if ( reordering_algo >= MAP ) {
      assert(map_file != "");
  }
  cout << "============ PARAMETERS ==========" << endl;
  cout << "threads: " << threads << " omp_get_max_threads: " << omp_get_max_threads() << endl;
  cout << "num_roots: " << num_roots << endl;
  cout << "map_file: " << map_file << endl;
  cout << "rand_gran: " << rand_gran << endl;
  cout << "reordering_algo: " << reordering_algo << endl;
  cout << "reordering_algo_str: " << ReorderingAlgoStr(reordering_algo) << endl;
  cout << "degree_used_for_reordering: " << degree_used_for_reordering << endl;
  cout << "degree_used_for_reordering_str: " << degree_used_for_reordering_str << endl;
  cout << "is_dense_write: " << isDenseWrite << endl;
  cout << "is_pagerank: " << isPageRank << endl;
#ifdef _SNIPER_
  cout << "sniper_min_threshold: " << sniper_min_threshold << endl;
#endif
  cout << "==================================" << endl;
#endif

  if (symmetric) {
      graph<symmetricVertex> G =
          readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
      pvector<uintE> new_ids(G.n, UINT_E_MAX);
      if ( reordering_algo != ORIGINAL ) {
#if 1 // Added by Priyank
          graph<symmetricVertex> newG;
          if ( reordering_algo > MAP ) {
              pvector<uintE> new_ids1(G.n, UINT_E_MAX);
              pvector<uintE> new_ids2(G.n, UINT_E_MAX);
              graph<symmetricVertex> newG1 = preprocessGraph<symmetricVertex>(G, symmetric, (degree_used_for_reordering == 0), new_ids1, false, false, MAP);
              newG = preprocessGraph<symmetricVertex>(newG1, symmetric, (degree_used_for_reordering == 0), new_ids2, false, false, (ReorderingAlgo)(reordering_algo - MAP));
              mergeTwoPreprocessingIndices(new_ids1, new_ids2, new_ids);
              newG1.del();
          } else {
              newG = preprocessGraph<symmetricVertex>(G, symmetric, (degree_used_for_reordering == 0), new_ids, false, false, reordering_algo);
          }
          G.del();
          Compute(newG,P,new_ids);
          for(int r=0;r<rounds;r++) {
              Compute(newG,P,new_ids);
          }
          newG.del();
#endif
      }
      else {
        Compute(G,P,new_ids);
        for(int r=0;r<rounds;r++) {
          Compute(G,P,new_ids);
        }
        G.del();
      }
    } else {
      graph<asymmetricVertex> G =
        readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
      pvector<uintE> new_ids(G.n, UINT_E_MAX);
      if ( reordering_algo != ORIGINAL ) {
#if 1 // Added by Priyank
          graph<asymmetricVertex> newG;
          if ( reordering_algo > MAP ) {
              pvector<uintE> new_ids1(G.n, UINT_E_MAX);
              pvector<uintE> new_ids2(G.n, UINT_E_MAX);
              graph<asymmetricVertex> newG1 = preprocessGraph<asymmetricVertex>(G, symmetric, (degree_used_for_reordering == 0), new_ids1, isPageRank, isDenseWrite, MAP);
              newG = preprocessGraph<asymmetricVertex>(newG1, symmetric, (degree_used_for_reordering == 0), new_ids2, isPageRank, isDenseWrite, (ReorderingAlgo)(reordering_algo - MAP));
              mergeTwoPreprocessingIndices(new_ids1, new_ids2, new_ids);
              newG1.del();
          } else {
              newG = preprocessGraph<asymmetricVertex>(G, symmetric, (degree_used_for_reordering == 0), new_ids, isPageRank, isDenseWrite, reordering_algo);
          }
        G.del();
        Compute(newG,P,new_ids);
        if(newG.transposed) newG.transpose();
        for(int r=0;r<rounds;r++) {
          Compute(newG,P,new_ids);
          if(newG.transposed) newG.transpose();
        }
        newG.del();
#endif
      }
      else {
        Compute(G,P,new_ids);
        if(G.transposed) G.transpose();
        for(int r=0;r<rounds;r++) {
          Compute(G,P,new_ids);
          if(G.transposed) G.transpose();
        }
        G.del();
      }
    }
#ifdef _SNIPER_
  if ( single_roi_nooutput != 0 && single_roi_output != 0 ) {
      cout << "ROI was nevery triggered for any iteration." << endl;
      cout << "Try setting sniper_min_threshold to a low number." << endl;
  }
#endif
}
#endif
