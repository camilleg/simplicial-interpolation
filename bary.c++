#include <algorithm>
#include <cmath>
#include <cstdio>

#include "bary.h"
#undef TESTING

static auto d = -1;

// Solve Ax=b, where A is d*d.
// Return false iff A is singular.
// Make a local copy of A.
// Crout's algorithm replaces that with its LU decomposition.
// Back-substitution solves for x, which is stored in b.
bool solveMatrix(const void* A, vertex& b)
{
  using dd = double[d][d];
  dd* src = (dd*)A;
  dd a;

  int i, j, k;
  for (i=0; i<d; ++i)
  for (j=0; j<d; ++j)
    a[i][j] = (*src)[i][j];

  double scaling[d];
  for (i=0; i<d; ++i)
    {
    double big = 0.;
    for (j=0; j<d; ++j)
      {
      const double t = fabs(a[i][j]);
      if (t > big)
        big = t;
      }
    if (big == 0.)
      return false; // a row of all zeros
    scaling[i] = 1. / big;
    }

  int permutation[d];
  int imax = -1;
  for (j=0; j<d; ++j)
    {
    for (i=0; i<j; ++i)
      {
      double t = a[i][j];
      for (k=0; k<i; ++k)
        t -= a[i][k] * a[k][j];
      a[i][j] = t;
      }
    double big = 0.;
    for (i=j; i<d; ++i)
      {
      double t = a[i][j];
      for (k=0; k<j; ++k)
	t -= a[i][k] * a[k][j];
      a[i][j] = t;
      t = scaling[i] * fabs(t);
      if (t >= big)
	{
	big = t;
	imax = i;
	}
      }
    if (j != imax)
      {
      // Swap rows.
      for (k=0; k<d; ++k)
	{
	const double t = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = t;
	}
      scaling[imax] = scaling[j];
      }
    permutation[j] = imax;
    if (a[j][j] == 0.)
      return false;
    if (j != d)
      {
      const double t = 1. / a[j][j];
      for (i=j+1; i<d; ++i)
	a[i][j] *= t;
      }
    }

  // Back substitution.
  int i1 = -1;
  for (i=0; i<d; ++i)
    {
    const int ip = permutation[i];
    double t = b[ip];
    b[ip] = b[i];
    if (i1 >= 0)
      for (int j=i1; j <= i-1; ++j)
        t -= a[i][j] * b[j];
    else if (t != 0.)
      i1 = i;
    b[i] = t;
    }
  for (i=d-1; i >= 0; --i)
    {
    double t = b[i];
    for (int j=i+1; j<d; ++j)
      t -= a[i][j] * b[j];
    b[i] = t / a[i][i];
    }
  return true;
}

// Calculate the determinant of the n*n matrix m, by recursing
// on "expansion by minors," i.e. discarding one row and column.
/*
  TODO:
  Convert this from recursion to iteration, by building a stack
  of n1 n1*n1 matrices.  Avoids function-call overhead.

  Would diagonalizing the matrix and multiplying its eigenvalues
  be any faster?  (Convert it to upper triangular form, avoiding
  numerical instabilities via Gauss-Jordan, then multiply the diagonal entries.)
  e.g., CLHEP 1.8.0.0 Matrix/Matrix.cc HepMatrix::dfact_matrix()
  That's O(n^3), same as what we have here.
*/
double determinant(const double* m, int n) {
  if (n < 1)
    return std::numeric_limits<double>::signaling_NaN();
  if (n == 1)
    return m[0];
  if (n == 2)
    return m[0] * m[3] - m[1] * m[2];
  /* if (n==3)
    return aei + cdh + bfg - ceg - bdi - afh; where abc is first row.  */

  const auto n1 = n-1;
  double m1[n1 * n1]; // Minor matrices.
  auto det = 0.0;
  for (auto j1=0; j1<n; ++j1) {      // For each row of m...
    // Compute m1.
    auto k=n, k1=0;
    for (auto i=1; i<n; ++i)        // Skip first column, i==0.
    for (auto j=0; j<n; ++j,++k)
      if (j != j1)                  // Skip j1'th row.
        m1[k1++] = m[k];            // k == n*i + j
    const auto t = determinant(m1, n1);

    // Recurse on m1.  Accumulate.
    det += m[j1] * (j1%2 ? -t : t);
  }
  return det;
}

simplexHint precomputeBary(const d_simplex& s, const vertex& centroid,
  const vector<vertex>& rgv, const vertex* raysCentroid, [[maybe_unused]] bool fRaySimplex)
{
  d = centroid.size();
  simplexHint h = {}; // h.s == nullptr, so h is (not yet) valid.
  if (d <= 1)
    return h;

#ifdef TESTING
  {
    d_simplex tmp(s);
    // s was sorted by sort_output().
    if (std::unique(tmp.begin(), tmp.end()) != tmp.end()) {
      dump_simplex("Internal error: duplicate vertices in simplex:", s);
      return h;
      }
  }
#endif

  // Precompute some properties of each facet.
  // ifacet'th facet is opposite ifacet'th vertex.
  resize(h.facetnormal, d, d+1);
  h.facetvertex.resize(d+1);
  h.facetvolume.resize(d+1);
  const auto fudge = pow(10.0, d*5-4); // Avoid astronomical volumes when d>3.
  for (int ifacet=0; ifacet<d+1; ++ifacet)
    {
    // Compute the inward-pointing unit-length normal to the facet.

    /*
    Given points v0, v1, ..., vn-1 (rgvFacet) in R^n,
      defining a facet of a (nondegenerate) n-simplex,
    find a vector w normal to the hyperplane that they span.
     
    Solution:
     
    Translate the hyperplane so v0 = 0.  That doesn't change the normal.
    (Computationally: subtract v0 from each of v1, ..., vn-1; then v0 := zero.
    Each row of A[d][d] is rgvFacet[i]-rgvFacet[0].)
      (Were it all zeros, A would be singular.)
     
    This reduces the problem to:
    find w simultaneously orthogonal to v1, ..., vn-1.
     
    i.e., w dot vi = 0, for i from 1 to n-1.
    n-1 equations in n unknowns.
    The other equation is simply something which {0} doesn't satisfy.
    Here, we arbitrarily choose "sum of w_i = 1."
     
    Normalize it to unit length.
    To orient it inwards:
    dot it with (centroid(s) - one_of_the_facet_points).
      If that's zero, the simplex is degenerate.
      If that's negative, negate it.
    */

    int i, j;

    // Compute the array of vertices of the facet by skipping
    // the ifacet'th vertex, the vertex opposite the facet.
    const vertex* rgvFacet[d];
#ifdef TESTING
    int schmoo[d];
#endif
    for (i=0,j=0; i<d+1; ++i)
      {
      const int ivertex = s[i];
#ifdef TESTING
      if (ivertex < 0 && !fRaySimplex)
        printf("internal error: inconsistency #1 in precomputeBary.\n");
      if (ivertex < 0 && i != d)
        printf("internal error: inconsistency #2: %d %d in precomputeBary.\n", i, d);
      if (fRaySimplex && i == d && ivertex != -1)
        printf("internal error: inconsistency #3: %d %d in precomputeBary.\n", i, d);
      schmoo[j] = ivertex;
#endif
      rgvFacet[j] = ivertex<0 ? raysCentroid : &rgv[ivertex];
      if (i != ifacet)
        if (++j >= d)
	  break;
      }

#ifdef TESTING
      for (j=0; j<d; ++j)
	printf("facet %d, rgvFacet[%d] = %d\n", ifacet, j, schmoo[j]);
      printf("\n");
#endif

    vertex& vNormal = h.facetnormal[ifacet];
    // vN = [1,0,...] and A's top row = [1,...] is the arbitrary equation
    // that vN's weights sum to unity.
    // If that happens to make A singular, we'll tweak A.
    vNormal[0] = 1.0;
    double a[d][d];
    for (j=0; j<d; ++j)
      a[0][j] = 1.0;
    for (i=1; i<d; ++i)
      for (j=0; j<d; ++j)
	a[i][j] = (*rgvFacet[i])[j] - (*rgvFacet[0])[j];

    if (d==2) {
      // Direct calculation.  (x,y) is v1-v0.
      // Agrees with solveMatrix().
      const auto x = (*rgvFacet[1])[0] - (*rgvFacet[0])[0];
      const auto y = (*rgvFacet[1])[1] - (*rgvFacet[0])[1];
      vNormal[0] = -y;
      vNormal[1] = x;
    } else {
      if (!solveMatrix(&a, vNormal)) {
	// Retry, with A tweaked to be nonsingular.  It *might* still
	// be singular, but let's wait and see if we need to be fussy.
	a[0][0] *= -1.0;
	std::fill(vNormal.begin()+1, vNormal.end(), 0.0);
	vNormal[0] = 1.0;
	if (!solveMatrix(&a, vNormal)) {
	  printf("precomputeBary(): matrix still singular.\n");
	  // for (i=0; i<d; ++i) { for (j=0; j<d; ++j) printf("%.1f\t", a[i][j]); printf("\n"); } printf("\n");
	  dump_simplex("possibly degenerate simplex:", s); // Degenerate faces?
	  return h;
	}
      }
    }

    // Scale vNormal to unit length.
    {
      auto t = 0.0;
      for (auto x: vNormal) t += sq(x);
      t = sqrt(t);
      for (auto& x: vNormal) x /= t;
    }

    // Compute vNormal dot (centroid - one_of_the_facet_points).
    {
      auto t = 0.0;
      for (i=0; i<d; ++i)
	t += vNormal[i] * (centroid[i] - (*rgvFacet[0])[i]);
      if (fabs(t) < 0.0001)
	printf("Warning: found degenerate simplex.\n");
      else if (t < 0.0)
	// Reverse normal to orient it inwards w.r.t. the simplex.
	for (auto& x: vNormal) x *= -1.0;
    }

#ifdef TESTING
    {
    // Verify that vNormal is orthogonal to each row of the (original) matrix.
    int i,j;

    for (int j=0; j<d; ++j)
      a[0][j] = 1.;
    for (int i=1; i<d; ++i)
      for (int j=0; j<d; ++j)
	a[i][j] = (*rgvFacet[i])[j] - (*rgvFacet[0])[j];

    for (i=1; i<d; ++i)
      {
      // Verify that a[i][] dot vNormal[] == 0.
      printf("Checking facet %d, row %d.\n", ifacet, i);
      double t = 0.;
      for (j=0; j<d; ++j)
        t += a[i][j] * vNormal[j];
      if (fabs(t) > .0001)
	{
        printf("Internal error in precomputeBary: vNormal isn't normal to facet %d's row %d:\n", ifacet, i);
	printf("  (%g, %g) dot (%g, %g) = %g, not zero.\n",
	    a[i][0], a[i][1], vNormal[0], vNormal[1], t);
	return h;
	}
      }
    }
#endif

    // Store a representative point on the facet's span.  We're just picking
    // a vertex of the facet, though the facet's centroid would increase
    // numerical stability if we want to get fancy.
    h.facetvertex[ifacet] = rgvFacet[0];

    // Optimization: if we cached/hashed these, we'd avoid computing each facet twice.

    // Compute the (d-1)-volume of each facet rgvFacet[0,...,d-1].
    auto& vol = h.facetvolume[ifacet];
    switch (d)
      {
    case 2:
	{
	// Fast: Euclidean distance from rgvFacet[0] to rgvFacet[1].
	const auto& a = *rgvFacet[0];
	const auto& b = *rgvFacet[1];
	vol = hypot(b[0]-a[0], b[1]-a[1]);
	break;
	}

    case 3:
	{
	// Fast: Area of triangle.
	// Vertices.
	const auto& a = *rgvFacet[0];
	const auto& b = *rgvFacet[1];
	const auto& c = *rgvFacet[2];
	// Side lengths.
	const auto la = sqrt(sq(b[0]-c[0]) + sq(b[1]-c[1]) + sq(b[2]-c[2]));
	const auto lb = sqrt(sq(c[0]-a[0]) + sq(c[1]-a[1]) + sq(c[2]-a[2]));
	const auto lc = sqrt(sq(b[0]-a[0]) + sq(b[1]-a[1]) + sq(b[2]-a[2]));
	// Heron's formula.
	const auto s = (la + lb + lc) * 0.5;
	vol = sqrt(s * (s-la) * (s-lb) * (s-lc));
	break;
	}

    default:
      // Use a generalization of Heron's formula, the Cayley-Menger determinant.
      // This also works for d=2 and d=3, but it's slow.
      double m[(d+1)*(d+1)];
      // Optimization:  m[] is symmetric in its main diagonal,
      // so we could compute only one half and reflect it into the other.
      int i,j;
      for (i=0; i<d+1; ++i)
        {
	m[i] = 1.;          // top row
	m[(d+1) * i] = 1.;  // left column
	}
      m[0] = 0.; // top left corner
      for (i=0; i<d; ++i)
      for (j=0; j<d; ++j)
        {
	// (i,j)-edge length squared
	double _ = 0.;
	if (i != j)
	  {
	  const auto& a = *rgvFacet[i];
	  const auto& b = *rgvFacet[j];
	  for (int k=0; k<d; ++k)
	    _ += sq(a[k] - b[k]);
	  }
	m[(i+1)*(d+1)+(j+1)] = _;
	}

      // Sloane's sequence A055546, (-1)^(d+1) / (2^d * (d!)^2).
      // More values, for up to d=100, are at https://oeis.org/A055546/b055546.txt.
      constexpr double scalars[16] = {
        -1, 2, -16, 288, -9216, 460800, -33177600, 3251404800, -416179814400,
        67421129932800, -13484225986560000, 3263182688747520000,
	-939796614359285760000.0, 317651255653438586880000.0,
        124519292216147926056960000.0, 56033681497266566725632000000.0 };
      const auto scalar = (d-1 < 16) ? scalars[d-1] : 1e40;
	// For huge d, approximation is okay, up to float overflow or underflow.
	// Only relative volumes matter because we'll normalize them anyways (to sum to 1).

      vol = sqrt(determinant(m, d+1) / scalar) / fudge;
      break;
      }
#ifdef TESTING
    if (vol < 0.0)
      printf("internal error in computing facet volumes.\n");
#endif
    }
  h.s = &s;
  return h;
}

// Fill w[0 to d+1] with q's barycentric coords w.r.t. s (i.e., w.r.t. h).
// Return true iff w[] is all nonnegative (i.e., q is in s.)
// If fRaySimplex, h must be in hiRay[], else hi[].
// If fRaySimplex, return true iff all w[!=d] are nonnegative.
bool computeBary(const simplexHint& h, const vertex& q, double* w, bool fRaySimplex)
{
  bool fPositive = true;
  int i;
  for (i=0; i<d+1; ++i)
    {
    // Dot product of (q - h.facetvertex[i]) with h.facetnormal[i].
    w[i] = 0.;
    for (int j=0; j<d; ++j)
      {
      const double dist = (q[j] - (*h.facetvertex[i])[j]) * h.facetnormal[i][j];
      w[i] += dist;
      }
    w[i] *= h.facetvolume[i];
    }

  {
    // Normalize w[i]'s so they sum to 1.
    auto sum = 0.0;
    for (i=0; i<d+1; ++i) sum += w[i];
    if (fabs(sum) < 1e-6) {
      // Typically for q with coords -1e20.  Fall back to "1/n" values, to avoid DBZ.
      printf("\tInternal workaround for zero-sum weights.\n");
      for (i=0; i<d+1; ++i)
	w[i] = copysign(1.0, w[i]);
      sum = 0.0;
      for (i=0; i<d+1; ++i) sum += w[i];
      if (fabs(sum) < 1e-6) {
	// Harsher fallback, to ensure nonzero sum.
	for (i=0; i<d+1; ++i) w[i] = 1.0;
	sum = d+1;
      }
    }
    for (i=0; i<d+1; ++i) {
      w[i] /= sum;
      if (w[i] < 0.0 && (!fRaySimplex || i != d))
	fPositive = false;
    }
  }

  if (fRaySimplex && fPositive && w[d] >= 0.)
    {
    printf("\tInternal inconsistency in raysimplex: no negative barycoord:\n\t");
    for (i=0; i<d+1; ++i) printf("%g ", w[i]);
    printf("\n");
    return false;
    }

  return fPositive;
}

/*
To compute the barycentric coordinates of a point q w.r.t. a simplex s:
Brute force:
    i'th coord := volume(s with i'th vertex replaced by q) / volume(s)
where all volumes are computed via the Cayley-Menger determinant
<http://mathworld.wolfram.com/Cayley-MengerDeterminant.html>.

Faster method inspired by H. Edelsbrunner:

q's barycentric coordinates are, by definition, the relative distances of q
from the facets of s.  Precompute the normals of all facets of s; then each
relative distance (each coordinate) is computed with only two scalar products.

Like this:
  abcde... are the facet's vertices.
  length := (q-a) dot normal, where normal has unit length
  and points "inwards into the simplex" from the facet.
  That way the barycentric coordinates
  will be negative if q lies outside the simplex, too.
  "a" could be any vertex of the facet, or any point on the facet's span.
*/
