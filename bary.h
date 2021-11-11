#include "si.h"

extern bool precomputeBary(const simplex& s, simplexHint& h, const vertex& centroid,
  const vertex* rgv, const vertex* raysCentroid, bool fRaySimplex = false);

extern bool computeBary(const simplexHint& h, const vertex& q, double* w,
  bool fRaySimplex = false);
