// Extend an R^d to R^e map from pointwise to continuous, via simplicial interpolation.

#pragma once
#include <array>
#include <iostream>
using std::cout;

constexpr auto d = 2;
constexpr auto e = 8;

#undef TESTING

// d-vertex and e-vertex.
// Most computation is done in d-space, so d is often implicit.
template <int n> using any_vertex = std::array<double, n>;
using   vertex = any_vertex<d>;
using e_vertex = any_vertex<e>;

template <class T> void dump_v(const char* prefix, const T& v) {
  cout << prefix;
  for (auto i: v) cout << i << " ";
  cout << "\n";
}

// Stores indices into the d-vertices qi[].
// -1 indicates qC, the common center point of the ray-simplices.
using d_simplex = std::array<int, d+1>;

inline void dump_simplex(const char* prefix, const d_simplex& s) {
  cout << prefix << "\n";
  extern vertex* qi;
  for (auto i: s) dump_v("\t", qi[i]);
}

// Inward-pointing normal vector of, and a point on, each facet of a simplex.
// Precomputing these speeds up computation of barycentric coordinates.
struct simplexHint {
  vertex facetnormal[d+1];
  const vertex* facetvertex[d+1];
  double facetvolume[d+1];
};

template <class T> T sq(const T& _) { return _*_; }
