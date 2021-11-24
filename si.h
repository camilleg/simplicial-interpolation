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
  constexpr int tSize[std::tuple_size<T>::value]{};
  cout << prefix;
  for (auto i=0u; i<std::size(tSize); ++i) cout << v[i] << " ";
  cout << "\n";
}

extern vertex* qi; // for d_simplex::dump()

// The x[]'s are indices into the vertices qi[];
// -1 indicates qC, the common center point of the ray-simplices.
// Stored in si[].
struct d_simplex
  {
  int x[d+1];
  int& operator[](int i) { return x[i]; }
  int  operator[](int i) const { return x[i]; }
  void dump(const char* sz = "s:") const
    { cout << sz << "\n"; for (auto i=0; i<d+1; ++i) dump_v("\t", qi[x[i]]); }
  };

// Inward-pointing normal vector of, and a point on, each facet of a simplex.
// Precomputing these speeds up computation of barycentric coordinates.
struct simplexHint {
  vertex facetnormal[d+1];
  const vertex* facetvertex[d+1];
  double facetvolume[d+1];
};

template <class T> T sq(const T& _) { return _*_; }
