#include "si.h"

bool precomputeBary(const simplex& s, simplexHint& h, const vertex& centroid, const vertex* rgv, const vertex* raysCentroid, bool fRaySimplex = false);

bool computeBary(const simplexHint& h, const vertex& q, double* w, bool fRaySimplex = false);
