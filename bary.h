// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

extern bool precomputeBary(const simplex& s, simplexHint& h, const vertex& centroid,
  const vertex* rgv, const vertex* raysCentroid, bool fRaySimplex = false);

extern bool computeBary(const simplexHint& h, const vertex& q, double* w,
  bool fRaySimplex = false);
