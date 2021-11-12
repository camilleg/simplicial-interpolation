#include <cmath>
#include <cstdio>

#include "si.h"
#include "det.h"

// Solve Ax=b, where A is d*d.  Returns false iff A is singular.
// Crout's algorithm replaces A with its LU decomposition.
// Then back-substitution solves for x, which is stored in b.
bool solveMatrix(double a[][d], vertex& b)
{
  int i, j, k;
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

bool precomputeBary(const simplex& s, simplexHint& h, const vertex& centroid,
  const vertex* rgv, const vertex* raysCentroid, [[maybe_unused]] bool fRaySimplex)
{
  if (d <= 1)
    return false;

#ifdef TESTING
  {
  for (int i=0; i<d+1; ++i)
    for (int j=i+1; j<d+1; ++j)
      if (s[i] == s[j])
        {
	printf("Internal error: simplex has duplicate vertices %d and %d\n", i, j);
	return false;
	}
  }
#endif

  // Precompute some properties of each facet.
  // ifacet'th facet is opposite ifacet'th vertex.
  for (int ifacet=0; ifacet<d+1; ++ifacet)
    {
    // Compute the inward-pointing unit-length normal to the facet.

    /*
    Given points v0, v1, ..., vn-1 in R^n,
      defining a facet of a (nondegenerate) n-simplex,
    Find a vector w normal to the hyperplane they span.
     
    Solution:
     
    Translate the plane so v0 = 0.  That won't change the normal.
    (Computationally: subtract v0 from each of v1, ..., vn-1; then v0 := zero.)
     
    This reduces the problem to:
    find w simultaneously orthogonal to v1, ..., vn-1.
     
    i.e., w dot vi = 0, for i from 1 to n-1.
    n-1 equations in n unknowns.
    The other equation is simply something which {0} doesn't satisfy.  e.g. sum of w_i = 1.
     
    Normalize it to unit length.
    To orient it inwards:
    dot it with (centroid(s) - one_of_the_facet_points).
      If that's zero, the simplex is degenerate.
      If that's negative, negate it.
    */

    int i, j;

    // Compute the array of vertices of the facet, not of the simplex.
    // By skipping the ifacet'th vertex, since that is the one opposite the facet.
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
      if (j >= d)
        printf("internal error in precomputeBary!\n");
      rgvFacet[j] = ivertex<0 ? raysCentroid : &rgv[ivertex];
      if (i != ifacet)
        if (++j >= d)
	  // Already finished;  quit now without overflowing.
	  break;
      }

#ifdef TESTING
      for (j=0; j<d; ++j)
	printf("facet %d, rgvFacet[%d] = %d\n", ifacet, j, schmoo[j]);
      printf("\n");
#endif

    vertex& vNormal = h.facetnormal[ifacet];
    for (i=0; i<d; ++i)
      vNormal[i] = i==0 ? 1. : 0.;

    double a[d][d];
    for (j=0; j<d; ++j)
      a[0][j] = 1.;
    for (i=1; i<d; ++i)
      for (j=0; j<d; ++j)
	{
	a[i][j] = (*rgvFacet[i])[j] - (*rgvFacet[0])[j];
	}

    if (!solveMatrix(a, vNormal))
      {
      printf("Internal error in precomputeBary().\n");
      // Matrix a was probably singular.  Because s had degenerate faces.
      s.dump("possibly degenerate simplex:");
      return false;
      }

    // Normalize vNormal to unit length (unfortunate pun).
    {
    double t = 0.;
    for (i=0; i<d; ++i)
      t += sq(vNormal[i]);
    t = sqrt(t);
    for (i=0; i<d; ++i)
      vNormal[i] /= t;
    }

    // Compute vNormal dot (centroid - one_of_the_facet_points).
    double t = 0.;
    for (i=0; i<d; ++i)
      t += vNormal[i] * (centroid[i] - (*rgvFacet[0])[i]);
    if (fabs(t) < .0001)
      printf("Warning: found degenerate simplex.\n");
    else if (t < 0.)
      // Reverse normal to orient it inwards w.r.t. the simplex.
      for (i=0; i<d; ++i)
	vNormal[i] *= -1.;

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
	return false;
	}
      }
    }
#endif

    // Store a representative point on the facet's span.  We're just picking
    // a vertex of the facet, though the facet's centroid would increase
    // numerical stability if we want to get fancy.
    h.facetvertex[ifacet] = rgvFacet[0];

    // Optimization: if we cached/hashed these, we'd avoid computing each facet twice.

    // Compute (d-1)-volume of each facet (defined by rgvFacet[0,...,d-1] ).
    switch (d)
      {
    case 2:
	{
	// Fast: Euclidean distance from rgvFacet[0] to rgvFacet[1].
	const auto& a = *rgvFacet[0];
	const auto& b = *rgvFacet[1];
	h.facetvolume[ifacet] = hypot(b[0]-a[0], b[1]-a[1]); // these values are OK, 8/29/02
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
	const double la = sqrt(sq(b[0]-c[0]) + sq(b[1]-c[1]) + sq(b[2]-c[2]));
	const double lb = sqrt(sq(c[0]-a[0]) + sq(c[1]-a[1]) + sq(c[2]-a[2]));
	const double lc = sqrt(sq(b[0]-a[0]) + sq(b[1]-a[1]) + sq(b[2]-a[2]));
	// Heron's formula.
	const double s = (la + lb + lc) * .5;
	h.facetvolume[ifacet] = sqrt(s * (s-la) * (s-lb) * (s-lc));
	break;
	}
    default:
      // Compute volume via a generalization of Heron's formula,
      // the Cayley-Menger determinant.

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

      const double scalars[14] = {
        -1.,2.,-16.,288.,-9216.,460800.,-33177600.,3251404800.,-416179814400.,
        67421129932800.,-13484225986560000.,3263182688747520000.,
	-939796614359285760000.,317651255653438586880000. };
      double scalar;
      if (d-1 <= 13)
        scalar = scalars[d-1];
      else
	{
        // Should really compute Sloane's sequence A055546: (-1)^(j+1) / (2^j * (j!)^2).
	// But approximation is okay (up to float overflow or underflow):
	// relative, not absolute, volumes is what matters
	// since we'll normalize them anyways (to sum to 1).
	scalar = 1e30; // vague approximation
	}

      h.facetvolume[ifacet] = sqrt(determinant(m, d+1) / scalar);
      break;
      }
#ifdef TESTING
    if (h.facetvolume[ifacet] < 0.)
      printf("internal error in computing facet volumes.\n");
#endif
    }
  return true;
}

// Fill w[0 to d+1] with q's barycentric coords w.r.t. s (i.e., w.r.t. h).
// Return true iff w[] is all nonnegative (i.e., q is in s.)
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
    double t = 0.;
    for (i=0; i<d+1; ++i)
      t += w[i];
    for (i=0; i<d+1; ++i)
      {
      w[i] /= t;
      if (w[i] < 0. && (!fRaySimplex || i!=d))
	fPositive = false;
      }
  }

  if (fRaySimplex && fPositive && w[d] >= 0.)
    {
    printf("internal inconsistency in raysimplex:  ALL positive barycoords.\n");
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
