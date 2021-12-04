#include "si.h"

simplexHint precomputeBary(const d_simplex& s, const vertex& centroid, const std::vector<vertex>& rgv, const vertex* raysCentroid, bool fRaySimplex);

bool computeBary(const simplexHint& h, const vertex& q, double* w, bool fRaySimplex = false);
