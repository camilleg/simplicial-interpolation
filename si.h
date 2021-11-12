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
using   vertex = std::array<double, d>;
using e_vertex = std::array<double, e>;
inline void dump_d(const char* sz, const   vertex& a) { cout << sz; for (int i=0; i<d; ++i) cout << a[i]; cout << "\n"; }
inline void dump_e(const char* sz, const e_vertex& a) { cout << sz; for (int i=0; i<e; ++i) cout << a[i]; cout << "\n"; }

extern vertex* qi; // for simplex::dump()

// A simplex's x[]'s are indices into the vertices qi[];
// -1 indicates qC, the common center point of the ray-simplices.
struct simplex
  {
  int x[d+1];
  int& operator[](int i) { return x[i]; }
  int  operator[](int i) const { return x[i]; }
  void dump(const char* sz = "s:") const
    { cout << sz << "\n"; for (auto i=0; i<d+1; ++i) dump_d("\t", qi[x[i]]); }
  };

// Inward-pointing normal vector of, and a point on, each facet of a simplex.
// Precomputing these speeds up computation of barycentric coordinates.
struct simplexHint
  {
  vertex facetnormal[d+1];
  const vertex* facetvertex[d+1];
  double facetvolume[d+1];
#ifdef UNUSED
  void dump(const char* sz = "hint:") const
    { cout << sz << "\n";
      for (int i=0; i<d+1; i++) dump_d("\tnormal: ", facetnormal[i]);
      for (int i=0; i<d+1; i++) dump_d("\tvertex: ", *facetvertex[i]);
      for (int i=0; i<d+1; i++) printf("\tvolume: %f\n", facetvolume[i]); }
#endif
  };

inline double sq(double _) { return _*_; }
