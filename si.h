// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

// Extending a pointwise map from R^d to R^e, to a continuous map
// via simplicial interpolation.

#pragma once

const int d = 2;
const int e = 8;

#undef TESTING

// (d-)vertex and e-vertex.
// Most computation is done in d-space, so d is often implicit.
class vertex
  {
 public:
  double x[d];
  double& operator[](int i) { return x[i]; }
  double  operator[](int i) const { return x[i]; }
  void dump(const char* sz = "v: ") const
    { printf("%s", sz); for (int i=0; i<d; i++) printf("%g ", x[i]); printf("\n"); }
  };
class e_vertex
  {
 public:
  double x[e];
  double& operator[](int i) { return x[i]; }
  double  operator[](int i) const { return x[i]; }
  void dump(const char* sz = "v: ") const
    { printf("%s", sz); for (int i=0; i<e; i++) printf("%g ", x[i]); printf("\n"); }
  };

extern vertex* qi; // for simplex::dump()

// Entries x[] of a simplex are indices into the array of vertices qi[].
// If -1, an entry indicates qC, the common center point of the ray-simplices.
class simplex
  {
 public:
  int x[d+1];
  int& operator[](int i) { return x[i]; }
  int  operator[](int i) const { return x[i]; }
  void dump(const char* sz = "s:") const
    { printf("%s\n", sz); for (int i=0; i<d+1; i++) qi[x[i]].dump("\t"); }
  void dumpi(const char* sz = "s:") const
    { printf("%s: ", sz); for (int i=0; i<d+1; i++) printf("%d ", x[i]); printf("\n"); }
  };

// Inward-pointing normal vector of, and a point on, each facet of a simplex.
// Precomputing these speeds up computation of barycentric coordinates.
class simplexHint
  {
 public:
  vertex facetnormal[d+1];
  const vertex* facetvertex[d+1];
  double facetvolume[d+1];
  void dump(const char* sz = "hint:") const
    { printf("%s\n", sz);
      for (int i=0; i<d+1; i++) facetnormal[i].dump("\tnormal: ");
      for (int i=0; i<d+1; i++) facetvertex[i]->dump("\tvertex: ");
      for (int i=0; i<d+1; i++) printf("\tvolume: %f\n", facetvolume[i]); }
  };

inline double sq(double _)
  { return _*_; }
