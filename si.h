// Extend an R^d to R^e map from pointwise to continuous, via simplicial interpolation.

#pragma once
#include <iostream>
#include <vector>
using std::cout;

#undef TESTING

// d-vertex and e-vertex.
using vertex = std::vector<double>;

template <class T> void dump_v(const char* prefix, const T& v) {
  cout << prefix;
  for (auto i: v) cout << i << " ";
  cout << "\n";
}

// Stores indices into the d-vertices qi[].
// -1 indicates qC, the common center point of the ray-simplices.
using d_simplex = std::vector<int>; // size d+1

extern void dump_simplex(const char* prefix, const d_simplex&);

// Inward-pointing normal vector of, and a point on, each facet of a simplex.
// Precomputing these speeds up computation of barycentric coordinates.
struct simplexHint {
  const d_simplex* s;
  std::vector<vertex> facetnormal;        // size d+1
  std::vector<const vertex*> facetvertex; // size d+1
  std::vector<double> facetvolume;        // size d+1
};

template <class T> T sq(const T& _) { return _*_; }

template <typename T> void resize(T& vec, int dim, int n) {
  vec.resize(n);
  for (auto& v: vec) v.resize(dim);
}
