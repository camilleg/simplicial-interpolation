// Extend an R^d to R^e map from pointwise to continuous, using simplicial interpolation.

#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
using std::cout;
using std::vector;

#undef TESTING

// Either a d-vertex or an e-vertex.
using vertex = vector<double>;

// Stores indices into the d-vertices qi[], so of size d+1.
// -1 indicates qC, the common center point of the ray-simplices.
using d_simplex = vector<int>;

// Inward-pointing normal vector of, and a point on, each facet of a simplex.
// All vectors have size d+1.
// Precomputing these speeds up computation of barycentric coordinates.
struct simplexHint {
  const d_simplex* s;
  vector<vertex> facetnormal;
  vector<const vertex*> facetvertex;
  vector<double> facetvolume;
};

// For testing and demos.
extern vertex qC; // Constructed common point of the ray-simplices.
extern vector<vertex> qi;
extern vertex pC; // What qC maps to.
extern vector<vertex> pi;
extern vector<d_simplex> si, siRay;

// Utilities.

void randomSites(vector<vertex>&, int dim, int n, double scale);

template <typename T> T sq(const T& _) {
  return _*_;
}

template <typename T> void resize(T& vec, int dim, int n) {
  vec.resize(n);
  for (auto& v: vec) v.resize(dim);
}

void dump_simplex(const char* prefix, const d_simplex&);

template <class T> void dump_v(const char* prefix, const T& v) {
  cout << prefix << std::setprecision(6);
  for (auto i: v) cout << i << " ";
  cout << "\n";
}

// API.

// Make n random points in R^e, and then corresponding ones in R^d.
enum class qi_kind { random, spaced, sammonsMapping, geneticAlgorithm };
bool init(int d, int e, int n, qi_kind kind);

// Map a point q in R^d to a point in R^e.
vertex eval(
  const vertex& q,
  bool* pfInside=nullptr,
  d_simplex* ps=nullptr,
  vertex* pcoords=nullptr,
  vertex* prFound=nullptr);

// Clean up.
void terminate();
