// Simplicial interpolation.

#include <algorithm>
#include <random>

#include "bary.h"
#include "edahiro.h"
#include "gacli.h"
#include "sammon.h"

#include "si.h"

vector<vertex> qi; // Points in R^d.
vector<vertex> pi; // Points in R^e.
vertex qC{0}; // Constructed common point of the ray-simplices.
vertex pC{0}; // What qC maps to.
vector<d_simplex> si, siRay; // Simplicial complex on qi[].
vector<simplexHint> hi, hiRay;

void dump_simplex(const char* prefix, const d_simplex& s) {
  cout << prefix << "\n";
  for (auto i: s) dump_v("\t", i<0 ? qC : qi[i]);
}

void randomSites(vector<vertex>& vec, int dim, int n, double scale) {
  resize(vec, dim, n);
  static std::default_random_engine rng;
  static std::uniform_real_distribution<double> range(0.0, 1.0);
  static bool fSeeded = false;
  if (!fSeeded) {
    fSeeded = true;
    rng.seed(std::random_device{}());
  }
  for (auto& v: vec) for (auto& x: v) x = scale * range(rng);
}

vertex add(const vertex& v, const vertex& w) {
  vertex r(v);
  for (auto i=0u; i < v.size(); ++i) r[i] += w[i];
  return r;
}

vertex subtract(const vertex& v, const vertex& w) {
  vertex r(v);
  for (auto i=0u; i < v.size(); ++i) r[i] -= w[i];
  return r;
}

double magnitude2(const vertex& v) {
  auto m = 0.0;
  for (auto& x: v) m += sq(x);
  return m;
}

double magnitude(const vertex& v) {
  return sqrt(magnitude2(v));
}

void rescale(vertex& v, double z) {
  for (auto& x: v) x *= z;
}

void setmag(vertex& v, double z) {
  const auto m = magnitude(v);
  if (m != 0.0)
    rescale(v, z / m);
}

// Scale the inputs to the hull algorithm, which uses exact integer arithmetic.
// Outside [1e2, 1e7], ch.c++ suffers degeneracies and overshoots.
// (Hull does this itself too, with mult_up.)
constexpr auto hullScale = 1e6;

void spacedSites(vector<vertex>& vec, int dim, int n) {
  randomSites(vec, dim, n, hullScale);
  // Iterate to repel points from each other.
  // Iterate fewer times for large n, to run faster.
  const auto iMax = 500u + 70000u / double(std::max(1, n));
  for (auto i=iMax; i>0; --i) {
    const auto relax = i / double(iMax);
    for (auto& v: vec) {
      vertex force(dim); // A direction, not a point.
      for (const auto& w: vec) {
	if (w == v)
	  continue;
	auto diff = subtract(v, w);
	const auto inv = relax * 2.5 * sq(hullScale) / magnitude2(diff);
	if (inv < 0.04) // Skip tiny forces.
	  continue;
	setmag(diff, inv);
	force = add(force, diff);
      }
      v = add(v, force);
    }
    // Get the bbox of the v's.
    constexpr auto big = std::numeric_limits<double>::max();
    vertex vMins(dim,  big);
    vertex vMaxs(dim, -big);
    for (const auto& v: vec) {
      for (auto i=0; i<dim; ++i) {
	vMins[i] = std::min(vMins[i], v[i]);
	vMaxs[i] = std::max(vMaxs[i], v[i]);
      }
    }
    // Per dimension, stretch to a bbox of size (0.0, hullScale).
    for (auto& v: vec) {
      for (auto i=0; i<dim; ++i) {
	v[i] -= vMins[i];
	v[i] *= hullScale / (vMaxs[i] - vMins[i]);
      }
    }
  }
}

// Sort each simplex's vertices, to compare multiple runs for testing,
// and to simplify other testing.
bool d_simplex_compare(const d_simplex& a, const d_simplex& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}
void sort_output(vector<d_simplex>& rgs, bool fRay) {
  for (auto& s: rgs)
    std::sort(s.begin(), s.end() - (fRay ? 1 : 0));
  std::sort(rgs.begin(), rgs.end(), d_simplex_compare);
}

bool init(int d, int e, int cPoint, qi_kind kind) {
  if (e < d) {
    printf("error: e (%d) < d (%d).\n", e, d);
    return false;
  }
  if (cPoint < d+1) {
    printf("error: not enough points for even one simplex: %d points < d+1 (%d).\n", cPoint, d+1);
    return false;
  }

  // Make output sites p_i.
  randomSites(pi, e, cPoint, 1.0);

  // Choose corresponding input sites q_i.
  switch (kind)
    {
  case qi_kind::random:
    // Uncorrelated with the p_i.  Uniformly distributed.
    randomSites(qi, d, cPoint, hullScale);
    break;
  case qi_kind::spaced:
    // Uncorrelated with the p_i, but roughly equidistant.
    spacedSites(qi, d, cPoint);
    break;
  case qi_kind::sammonsMapping:
    resize(qi, d, cPoint);
    computeSammon(qi, pi, hullScale);
    break;
  case qi_kind::geneticAlgorithm:
    // Strongly correlated with the p_i.
    // If tmp were on the stack, if e*cPoint > 200000 or so, the stack might overflow.
    auto tmp = new double[e*cPoint];
    for (auto i=0; i<cPoint; ++i)
    for (auto j=0; j<e; ++j)
      tmp[i*e + j] = pi[i][j];
    const auto m = GADistanceMatrix(cPoint, e, d, tmp);
    delete [] tmp;
    resize(qi, d, cPoint);
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<d; ++j)
      qi[i][j] = hullScale * double(m->rgl[i*d + j]) / sHuge;
    free(m);
    break;
    }
  // randomSites() and computeSammon() make no errors from bogus input,
  // but delaunay_tri() reports degeneracies and Edahiro_Init() fails gracefully.
  // GADistanceMatrix() handles errors (from bogus input) simply by crashing.

  // Store a triangulation of the qi's in si.
  // callhull.c++ avoids tightly coupling hull.h to this file,
  // coupling Ken Clarkson's hull code to Camille's simplicial interpolation code.
  extern bool delaunay_tri(vector<d_simplex>&, vector<d_simplex>&, int, int);
  if (!delaunay_tri(si, siRay, d, qi.size()) || si.empty()) {
    printf("error: made no simplices.\n");
    // Maybe d or cPoint is too small.
    return false;
  }
  sort_output(si, false);
  sort_output(siRay, true);
#ifdef DUMP_SIMPLICES
  for (auto s: si) {
    for (auto i: s) printf("%2d ", i);
    printf("\n");
  }
  for (auto s: siRay) {
    for (auto i: s) printf("%2d ", i);
    printf("\n");
  }
#endif

  // printf("read %lu true simplices, %lu ray-simplices.\n", si.size(), siRay.size());

  // Output is a sequence of lines.  Each line lists the vertex indices of a d-simplex.
  // An index of -1 indicates the point at infinity;  we'll use that for
  // building ray-simplices.

  if (d == 2 && !Edahiro_Init(qi, si)) {
    printf("error: Edahiro failed to init.\n");
    return false;
  }

  // Precompute some things to speed up eval().

  // Calculate qC, the centroid of the bounding box of all simplices.
  qC.resize(d);
  for (auto i=0; i<d; ++i)
    {
    auto zMin = std::numeric_limits<double>::max();
    auto zMax = -zMin;
    for (const auto& v: qi)
      {
      const double z = v[i];
      if (z < zMin)
        zMin = z;
      if (z > zMax)
        zMax = z;
      }
    qC[i] = (zMin + zMax) * .5;
    }

  hi.clear();
  for (const auto& s: si) {
    // Accumulate into vC the centroid of s.  Pass that to precomputeBary().
    vertex vC(d);
    for (auto j=0; j<d; ++j) {
      for (auto i: s) vC[j] += qi[i][j]; // k<0 is possible only for siRay, not for si.
      vC[j] /= d + 1.0;
    }
    const auto h = precomputeBary(s, vC, qi, &qC, false);
    if (!h.s) {
      printf("error: precomputeBary failed.\n");
      return false;
    }
    hi.push_back(h);
  }

  hiRay.clear();
  for (const auto& s: siRay) {
    // Accumulate into vC the centroid of s.  Pass that to precomputeBary().
    vertex vC(d);
    for (auto j=0; j<d; ++j) {
      for (auto i: s) vC[j] += (i < 0 ? qC : qi[i])[j];
      vC[j] /= d + 1.0;
    }
    const auto h = precomputeBary(s, vC, qi, &qC, true);
    if (!h.s) {
      printf("error: ray-simplex precomputeBary failed.\n");
      return false;
    }
    hiRay.push_back(h);
  }

  pC = eval(qC);
  return true;
}

void terminate()
{
  pi.clear();
  qi.clear();
  si.clear();
  siRay.clear();
  hi.clear();
  hiRay.clear();
}

// Return the simplex that contains q.
// On error, return an arbitrary simplex.
// Into w[0...d+1], stuff q's barycentric coordinates w.r.t. that simplex.
const d_simplex& findSimplex(const vertex& q, double* w, bool* pfInside=nullptr) {
  if (pfInside) *pfInside = true;
  auto failover = false;
  if (q.size() == 2) {
    // Edahiro's algorithm.  Fast.
    const int i = Edahiro_RegionFromPoint(q[0], q[1]);
    if (i >= 0) {
      if (i >= int(si.size())) {
	printf("\tinternal error: edahiro returned out-of-range simplex %d > %lu\n", i, si.size());
	failover = true; goto Lfailover;
      }
      if (!computeBary(hi[i], q, w)) {
	printf("\tinternal error: edahiro returned wrong simplex %d from (%.1f, %.1f).\n", i, q[0], q[1]);
	failover = true; goto Lfailover;
      }
      return si[i];
    }
  } else {
    // Brute force.
    // Compute q's bary-coords w.r.t. each simplex in si[].
    // If one has coordinates all nonnegative, return that one.
Lfailover:
    for (const auto& h: hi)
      if (computeBary(h, q, w)) {
	if (failover) {
	  const auto it = std::find(si.begin(), si.end(), *h.s);
	  printf("\tedahiro should have returned %ld\n", std::distance(si.begin(), it));
	}
	return *h.s;
      }
  }
  // q wasn't in any simplex, so look in the ray-simplices.
  if (pfInside) *pfInside = false;
  for (const auto& h: hiRay)
    if (computeBary(h, q, w, true))
      return *h.s;
  // This should be impossible, because the ray-simplices partition R^d.
  // Arbitrarily return the first ray-simplex.
  printf("\tinternal error in findSimplex\n");
  (void)computeBary(hi[0], q, w);
  return siRay[0];
}

// Map a d-vertex to an e-vertex.
// These optional args are ugly.  A better design may appear after a few more demos have been written.
vertex eval(const vertex& q, bool* pfInside, d_simplex* ps, vertex* pcoords, [[maybe_unused]] vertex* prFound) {
  // Find which simplex s contains q.
  const auto d = q.size();
  double w[d+1]; // q's coordinates w_j with respect to s.
  const d_simplex& s = findSimplex(q, w, pfInside);

#ifdef TESTING
  // Verify that q == the point whose barycoords are w[] wrt s.
  vertex r(d);
  for (auto j=0u; j<d; ++j) {
    for (auto i=0u; i<d+1; ++i)
      r[j] += w[i] * (s[i] < 0 ? qC : qi[s[i]])[j];
  }
  // The reconstructed point is r.  How far is it from q?
  // (How accurate were the barycoords w[]?)
  auto dist = 0.0;
  for (auto j=0u; j<d; ++j)
    dist += sq(r[j] - q[j]);
  dist = sqrt(dist);
  if (dist > 1e-8)
    printf("warning: reconstruction error = %g\n\n\n", dist);
  if (prFound) *prFound = r;
#endif

  if (ps) *ps = s;
  if (pcoords) std::copy(w, w + d+1, pcoords->begin()); // For discs.

  // Sum with weights w[] and vertices pi[s[]].
  const auto e = pi.front().size();
  vertex p(e);
  pC = p; // Don't crash in "pC = eval(qC);"
  for (auto j=0u; j<e; ++j)
    for (auto i=0u; i<d+1; ++i)
      p[j] += w[i] * (s[i] < 0 ? pC : pi[s[i]])[j];
  return p;
}
