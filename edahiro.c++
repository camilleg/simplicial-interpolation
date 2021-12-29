// Constant-time lookup of which triangle contains a query point.
//
// Bucketing algorithm.  C++ translation of [Edahiro84]'s FORTRAN code.
// Names of identifiers are generally identical to the original FORTRAN.
//
// Vertices are numbered 1 to nv.
// Edges are numbered 1 to ne.
// Regions are numbered 1 to nr, and 0 for the unique unbounded region.

#include <cmath>
#include <cstdio>
#include <cstring> // memset

#include "edahiro.h"

struct Vertex {
  double x, y;
  int ies; // index into rgedge, an edge incident to this vertex
};

static constexpr auto cvertexMax = 8200;
static Vertex rgvertex[cvertexMax+1];
static auto nv = 1;

struct Edge {
  int ihead, itail; // ivertex's
  int lface, rface;
  int pnext, mnext; // iedges's: ccw-neighbor edges incident to ihead, itail
  double angle;     // "slope" [0,2PI) of forward edge
  double negangle;  // "slope" [0,2PI) of backward edge
  static constexpr auto NaN = std::numeric_limits<double>::signaling_NaN();

  Edge(): ihead(-1), itail(-1), lface(-1), rface(-1), pnext(-1), mnext(-1), angle(NaN), negangle(NaN) {}
};

struct Region {
  int i1, i2, i3; // ipt's (our regions are only triangular).
};
static constexpr auto cregionMax = 1000; // 350*250?
static Region rgregion[cregionMax+2];
static auto nr = 1;

static constexpr auto cedgeMax = 12300;
static Edge rgedge[cedgeMax+1];
static auto ne = 1;

inline double  vx(const int i) { return rgvertex[i].x; }
inline double  vy(const int i) { return rgvertex[i].y; }
inline int ies(const int i) { return rgvertex[i].ies; }

inline int ihead(const int i) { return rgedge[i].ihead; }
inline int itail(const int i) { return rgedge[i].itail; }
inline int lface(const int i) { return rgedge[i].lface; }
inline int rface(const int i) { return rgedge[i].rface; }
inline int pnext(const int i) { return rgedge[i].pnext; }
inline int mnext(const int i) { return rgedge[i].mnext; }
inline double angle(const int i) { return rgedge[i].angle; }
inline double negangle(const int i) { return rgedge[i].negangle; }

// Data structures for preprocessing phase.
static int ihor[351][251];
static int hor[30001][3];
static int iver[351][251];
static double ver[3001];
static int pver[30001][3];
static int inod[351][251];
static int nod[8201];
static int face[351][251];
static int edge[12301][2]; // Doubly linked list of edges.  [0] is prev, [1] is next.

static auto xmin = 0.0;
static auto ymin = 0.0;
static auto xminOrig = 0.0;
static auto yminOrig = 0.0;
static auto xd = 0.0;
static auto yd = 0.0;
static auto xdOrig = 0.0;
static auto ydOrig = 0.0;

static auto n1 = 0; // number of buckets
static auto n2 = 0;

inline bool Preprocess()
{
  int i;
  int nheh = 0;
  int nhev = 0;

  // Compute widths of the bucketed space, xd and yd.
  double xmax = vx(1);
  double ymax = vy(1);
  xmin = xmax;
  ymin = ymax;
  // Our caller Edahiro_Init verified that nv >= 3.
  for (i=2; i<=nv; i++)
    {
    if (vx(i) > xmax) xmax = vx(i);
    if (vx(i) < xmin) xmin = vx(i);
    if (vy(i) > ymax) ymax = vy(i);
    if (vy(i) < ymin) ymin = vy(i);
    }
  xdOrig = xd = xmax - xmin;
  ydOrig = yd = ymax - ymin;
  if (xd < 1.0 || yd < 1.0) {
    fprintf(stderr, "Preprocess failed: inputs too similar.\n");
    return false;
  }
  xminOrig = xmin;
  yminOrig = ymin;

  // Compute Sx and Sy (vnx and vny, respectively).
  double vnx = 0.;
  double vny = 0.;
  for (i=1; i<=ne; i++)
    {
    const int i1 = ihead(i);
    const int i2 = itail(i);
    vnx += fabs(vx(i1) - vx(i2));
    vny += vy(i1) - vy(i2); // No abs() because edges point upwards.
    }

  // Cover a frame over the given graph and partition it into buckets.

  vnx /= xd;
  vny /= yd;
  //printf("	xd = %.1f    yd = %.1f\n", xd, yd);
  //printf("	Sx = %f    Sy = %f\n", vnx, vny);

  const double vrate = 2.; // "alpha"  (from p.22 of preprint)
  const double vper = 1.5; // "beta"   (from p.22 of preprint)
  const int nsum = int(vrate * vrate * nv/2.0);
  const double vrate1 = vrate * sqrt(nv / vper);
  n1 = int(vrate1 * sqrt(vny / vnx) * vper + 1.);
  n2 = 0;
  if (n1 <= 2)
    {
    n1 = 2;
    n2 = nsum;
    }
  else
    {
    n2 = int(vrate1 * sqrt(vnx / vny) + 1.);
    if (n2 <= 2)
      {
      n2 = 2;
      n1 = nsum;
      }
    }

  //printf("\tBuckets: %d by %d.\n", n1, n2);
  xd /= (n1 * 4.);
  yd /= (n2 * 4.);
  xmin -= xd;
  ymin -= yd;
  xd *= 4. + 2./n1;
  yd *= 4. + 2./n2;

  // Zero the arrays.
  for (int iy=1; iy<=n2; iy++) {
    face[1][iy] = 0;
    for (int ix=1; ix<=n1; ix++) {
      ihor[ix][iy] = 0;
      iver[ix][iy] = 0;
      inod[ix][iy] = 0;
    }
  }

  // Transform coordinate system.

  for (i=1; i<=nv; i++)
    {
    rgvertex[i].x = (vx(i) - xmin) / xd + 1.;
    rgvertex[i].y = (vy(i) - ymin) / yd + 1.;
    const int ix = vx(i);
    const int iy = vy(i);

    // Construct Nij
    int i0 = inod[ix][iy];
    if (i0 == 0 || vy(i0) <= vy(i))
      {
      nod[i] = i0;
      inod[ix][iy] = i;
      continue;
      }
L101:
    const int i1 = nod[i0];
    if (i1 == 0 || vy(i1) <= vy(i))
      {
      nod[i0] = i;
      nod[i] = i1;
      continue;
      }

    i0 = i1;
    goto L101;
    }

  // Construct Vij and Hij.
  // yyd is the gradient of the edge, in the bucket coordinate system.

  for (i=1; i<=ne; i++)
    {
    edge[i][0] = 0;
    edge[i][1] = 0;
    int i1 = itail(i);
    int i2 = ihead(i);
    const int j1 = int(vx(i1));
    const int j2 = int(vx(i2));
    int iy = -1;
    int iy1;
    int j;
    double yj1;
    double yyd;

    if (j1 == j2)
      {
      // Edge crosses no vertical partition line.

      j = j1 + 1;
      iy  = int(vy(i1));
      iy1 = int(vy(i2));
      }
    else if (j1 < j2)
      {
      // Edge is directed upper-left.

      yyd = (vy(i1) - vy(i2)) / (vx(i1) - vx(i2));
      iy1 = int(vy(i1));
      yj1 = vy(i1) + yyd * (double(j1) - vx(i1));
      for (j = j1+1; j<=j2; j++)
	{
	yj1 += yyd;
	++nhev;
	pver[nhev][1] = i;
	ver[nhev] = yj1;
	iy = int(yj1);
	if (iy1 != iy)
	  {
	  for (int iyy = iy1; iyy <= iy-1; iyy++)
	    {
	    ++nheh;
	    hor[nheh][1] = i;
	    hor[nheh][2] = ihor[j-1][iyy];
	    ihor[j-1][iyy] = nheh;
	    }
	  }
	int i0 = iver[j][iy];
	if (i0 == 0 || ver[i0] <= yj1)
	  {
	  pver[nhev][2] = i0;
	  iver[j][iy] = nhev;
	  }
	else
	  {
	  for (;;)
	    {
	    i1 = pver[i0][2];
	    if (i1 == 0 || ver[i1] <= yj1)
	      break;
	    i0 = i1;
	    }
	  pver[i0][2] = nhev;
	  pver[nhev][2] = i1;
	  }
	iy1 = iy;
	}
      
      iy1 = int(vy(i2));
      }
    else // j1 > j2
      {
      // Edge is directed upper-right.

      yyd = (vy(i1) - vy(i2)) / (vx(i1) - vx(i2));
      iy1 = int(vy(i2));
      yj1 = vy(i2) + yyd * (double(j2) - vx(i2));
      for (j = j2+1; j <= j1; j++)
	{
	yj1 += yyd;
	++nhev;
	pver[nhev][1] = i;
	ver[nhev] = yj1;
	iy = int(yj1);
	if (iy1 != iy)
	  {
	  for (int iyy = iy; iyy <= iy1-1; iyy++)
	    {
	    ++nheh;
	    hor[nheh][1] = i;
	    hor[nheh][2] = ihor[j-1][iyy];
	    ihor[j-1][iyy] = nheh;
	    }
	  }
	int i0 = iver[j][iy];
	if (i0 == 0 || ver[i0] <= yj1)
	  {
	  pver[nhev][2] = i0;
	  iver[j][iy] = nhev;
	  }
	else
	  {
	  for (;;)
	    {
	    i2 = pver[i0][2];
	    if (i2 == 0 || ver[i2] <= yj1)
	      break;
	    i0 = i2;
	    }
	  pver[i0][2] = nhev;
	  pver[nhev][2] = i2;
	  }
	iy1 = iy;
	}
      
      iy = int(vy(i1));
      }

    if (iy != iy1)
      for (int iyy = iy; iyy <= iy1-1; iyy++)
	{
	hor[++nheh][1] = i;
	hor[  nheh][2] = ihor[j-1][iyy];
	ihor[j-1][iyy] = nheh;
	}
    }

  // Construct Fij.

  if (n1 != 1 && n2 != 1)
    for (int ix=2; ix<=n1; ix++)
      {
      face[ix][n2] = 0;
      int sface = 0;
      for (int iy=1; iy<=n2-1; iy++)
	{
	const int n0 = iver[ix][iy];
	if (n0 != 0)
	  {
	  const int e0 = pver[n0][1];
	  sface = (vx(itail(e0)) > vx(ihead(e0))) ? rface(e0) : lface(e0);
	  }
	face[ix][iy] = sface;
	}
      }

//printf("\n\t____ face[][]: ____\n\t");
//for (int iy=1; iy<=n2; iy++)
//  {
//  for (int ix=1; ix<=n1; ix++)
//    printf("%2d ", face[ix][iy]);
//  printf("\n\t");
//  }
//printf("\n");
  return true;
}

int DoQuery(double x, double y)
{
  // Transform (x,y) coordinates into bucket coordinates.
  //
  // The query point lies in Bij.  ipx is i, ipy is j.

  x = (x-xmin)/xd + 1.;
  const int ipx = int(x);
  y = (y-ymin)/yd + 1.;
  const int ipy = int(y);

  // If point is out of range, it must be the unbounded region.
  if (ipx<0 || ipy<0 || ipx>n1 || ipy>n2)
    return 0;

  // STEP ONE
  //
  // Doubly linked list "edge" represents S.  "top" is its header.

  int top = ne + 1;
  int i = ihor[ipx][ipy];
  while (i != 0)
    {
    const int i1 = hor[i][1];
    edge[i1][0] = top;
    edge[top][1] = i1;
    top = i1;
    i = hor[i][2];
    }
  int left = face[ipx][ipy];

  // Compute S = symmetric difference of S and Vij(Q).

  int i1 = iver[ipx][ipy];
  if (i1 != 0 && ver[i1] > y)
    {
    do
      {
      i = pver[i1][1];
      if (edge[i][0] == 0)
	{
	edge[i][0] = top;
	edge[top][1] = i;
	top = i;
	}
      else if (edge[i][1] == 0)
	{
	top = edge[i][0];
	edge[top][1] = 0;
	edge[i][0] = 0;
	}
      else
	{
	edge[edge[i][1]][0] = edge[i][0];
	edge[edge[i][0]][1] = edge[i][1];
	edge[i][0] = 0;
	edge[i][1] = 0;
	}
      i1 = pver[i1][2];
      }
    while (i1 != 0 && ver[i1] > y);

    // Assign f.

    left = vx(itail(i)) < vx(ihead(i)) ? rface(i) : lface(i);
    }

  // Compute S = symmetric difference of S and Nij(Q).

  int i0 = inod[ipx][ipy];
  while (i0 != 0 && vy(i0) > y)
    {
    const int i2 = i1 = ies(i0);
    do
      {
      i = i1>=0 ? i1 : -i1;
      if (edge[i][0] == 0)
	{
	edge[i][0] = top;
	edge[top][1] = i;
	top = i;
	}
      else if (edge[i][1] == 0)
	{
	top = edge[i][0];
	edge[top][1] = 0;
	edge[i][0] = 0;
	}
      else
	{
	edge[edge[i][1]][0] = edge[i][0];
	edge[edge[i][0]][1] = edge[i][1];
	edge[i][0] = 0;
	edge[i][1] = 0;
	}
      i1 = i1 < 0 ? mnext(i) : pnext(i);
      }
    while (i1 != i2);
    i0 = nod[i0];
    }

  // Compute S = symmetric difference of S and Vi+1,j(Q).

  i1 = iver[ipx+1][ipy];
  while (i1 != 0 && ver[i1] > y)
    {
    i = pver[i1][1];
    if (edge[i][0] == 0)
      {
      edge[i][0] = top;
      edge[top][1] = i;
      top = i;
      }
    else if (edge[i][1] == 0)
      {
      top = edge[i][0];
      edge[top][1] = 0;
      edge[i][0] = 0;
      }
    else
      {
      edge[edge[i][1]][0] = edge[i][0];
      edge[edge[i][0]][1] = edge[i][1];
      edge[i][0] = 0;
      edge[i][1] = 0;
      }
    i1 = pver[i1][2];
    }

  // STEP TWO

  i0 = 0;
  i = top;
  double x1 = ipx;
  while (i <= ne)
    {
    int i1 = ihead(i);
    const int i2 = itail(i);
    const double wx = vx(i1);
    const double wy = vy(i1);
    const double x2 = (y-wy) * (vx(i2)-wx) / (vy(i2)-wy) + wx;
    if (x1 < x2 && x > x2)
      {
      x1 = x2;
      i0 = i;
      }
    i1 = edge[i][0];
    edge[i][0] = 0;
    edge[i][1] = 0;
    i = i1;
    }

  // Compute the answer and return it.

  if (i0 == 0)
    return left;

  if (vy(itail(i0)) == y)
    {
    // The query point lies on the boundary of the slab.

    // Avoid infinite loops in some degenerate cases (roundoff
    // error), by arbitrarily choosing one of the iterated edges.
    int cloop = 0;

    int i1 = ihead(i0);
    double x0 = vx(i1) - x1;
    double y1 = vy(i1) - y;
    while ((i0 = pnext(i0)) >= 0)
      {
      i1 = ihead(i0);
      const double x2 = vx(i1) - x1;
      const double y2 = vy(i1) - y;
      if (x2 * y1 > x0 * y2 || ++cloop > ne+2)
	{
	if (cloop > ne+2)
	  printf("Edahiro warning: avoided infinite loop!\n");
	return rface(i0);
	}

      x0 = x2;
      y1 = y2;
      }

    while ((i0 = mnext(-i0)) < 0)
      ;
    }

  return rface(i0);
}

// Adjust angle to lie in the range [0, 2PI).
inline double Normalize(const double a)
  { return a<0. ? a+2.*M_PI : a>=2.*M_PI ? a-2.*M_PI : a; }

inline double Angle(const double dy, const double dx)
  { return Normalize(atan2(dy, -dx)); }

inline int Tweak(const int ir)
  { return ir >= 1 ? ir : 0; }
  // 0 is index of (implicit) infinite unbounded region

// Zero-based, not 1-based.
// Permuting of points so they are sorted by y-coordinate.

bool Edahiro_Init(const vector<vertex>& qi, const vector<d_simplex>& si)
{
  // Read the output of Hull algorithm.
  // nv, x's and y's, ntri, itri's.
  // Stuff nv and rgvertex, nr and rgregion.

  nv = qi.size();
  if (nv < 3) {
    fprintf(stderr, "error: Edahiro_Init() needs at least 3 vertices, not %d.\n", nv);
    return false;
  }
  nr = si.size();
  if (nr < 1) {
    fprintf(stderr, "error: Edahiro_Init() needs at least 1 region, not %d.\n", nr);
    return false;
  }

  for (auto iv=1; iv<=nv; ++iv) {
    rgvertex[iv].x = qi[iv-1][0];
    rgvertex[iv].y = qi[iv-1][1];
  }
  for (auto ir=1; ir<=nr; ++ir) {
    // Change indices from zero-based to one-based (FORTRAN-ism).
    rgregion[ir].i1 = si[ir-1][0] + 1;
    rgregion[ir].i2 = si[ir-1][1] + 1;
    rgregion[ir].i3 = si[ir-1][2] + 1;
    //if (ir<5) printf("new raw region: %d %d %d\n", rgregion[ir].i1, rgregion[ir].i2, rgregion[ir].i3);
  }

  // Sort rgvertex[] by increasing y-coordinate.
  // (Bubble-sort for now.)
  // Adjust rgregion[].* to match.
  // For each swap in rgvertex, swap the rgregion[].i_'s with those values too.
  for (auto iv=nv-2; iv>1; --iv) {
    for (auto jv=1; jv<iv; ++jv) {
      Vertex& r0 = rgvertex[jv];
      Vertex& r1 = rgvertex[jv+1];
      if (r0.y > r1.y) {
	std::swap(r0.x, r1.x);
	std::swap(r0.y, r1.y);
	// Swap all rgregion[].i_'s with the values jv and jv+1.
	for (auto ir=1; ir<=nr; ++ir) {
	  Region& r = rgregion[ir];
	  if (r.i1 == jv) { r.i1 = jv+1; } else if (r.i1 == jv+1) { r.i1 = jv; }
	  if (r.i2 == jv) { r.i2 = jv+1; } else if (r.i2 == jv+1) { r.i2 = jv; }
	  if (r.i3 == jv) { r.i3 = jv+1; } else if (r.i3 == jv+1) { r.i3 = jv; }
	}
      }
    }
  }

  //for (auto ir=1; ir<=nr; ++ir)
  //  printf("new region: %d %d %d\n", rgregion[ir].i1, rgregion[ir].i2, rgregion[ir].i3);

  {
    int rgLeftRegion [nv+1][nv+1];
    int rgRightRegion[nv+1][nv+1];
    memset(rgLeftRegion , -1, sizeof rgLeftRegion);
    memset(rgRightRegion, -1, sizeof rgRightRegion);

    // Stuff all the other fields.
    for (auto ir=1; ir<=nr; ++ir) {
      Region& r = rgregion[ir];

      // Sort r's vertices counterclockwise.
      // Because r is a triangle, just check one of its angles.
      // If that exceeds 180 degrees, go around it the other way.
      //if (ir<5) fprintf(stderr, "region %d %d %d\n", ir, r.i1, r.i2);
      const auto angleAB = Angle(vy(r.i2) - vy(r.i1), vx(r.i2) - vx(r.i1));
      const auto angleAC = Angle(vy(r.i3) - vy(r.i1), vx(r.i3) - vx(r.i1));
      if (Normalize(angleAC - angleAB) > M_PI)
	std::swap(r.i2, r.i3);

      // Record the regions to the left of these edges,
      // and to the right of the reversed edges.
      rgLeftRegion [r.i1][r.i2] = ir;
      rgLeftRegion [r.i2][r.i3] = ir;
      rgLeftRegion [r.i3][r.i1] = ir;
      rgRightRegion[r.i2][r.i1] = ir;
      rgRightRegion[r.i3][r.i2] = ir;
      rgRightRegion[r.i1][r.i3] = ir;

      //printf("leftrgn(%d) of %d %d %d\n", ir, r.i1, r.i2, r.i3);
      //printf("rigtrgn(%d) of %d %d %d\n", ir, r.i1, r.i3, r.i2);
    }

    // Accumulate the edge list.
    // We can hit only the lower triangle of rgLeftRegion,rgRightRegion
    // (jv<iv, not jv<=nv) because the ir-loop above hits each edge
    // in both directions.

    ne = 1;  // One-based FORTRANism.
    for (auto iv=1; iv<=nv; ++iv)
    for (auto jv=1; jv< iv; ++jv) {
      if (rgLeftRegion[iv][jv] <= 0 && rgRightRegion[iv][jv] <= 0)
	// No edge joins the iv'th and jv'th vertices.
	continue;

      Edge& e = rgedge[ne++];
      const bool edgeForwards = (vy(iv) > vy(jv) || (vy(iv)==vy(jv) && vx(iv)<vx(jv)));
      e.ihead = edgeForwards ? iv : jv;
      e.itail = edgeForwards ? jv : iv;
      e.lface = Tweak(rgLeftRegion [e.ihead][e.itail]);
      e.rface = Tweak(rgRightRegion[e.ihead][e.itail]);
      e.angle    = Angle(vy(e.itail) - vy(e.ihead), vx(e.itail) - vx(e.ihead));
      e.negangle = Angle(vy(e.ihead) - vy(e.itail), vx(e.ihead) - vx(e.itail));
    }
    --ne;
    // For all e above, vy(e.ihead) >= vy(e.itail);  if ==, vx(head) < vx(tail).
    // This is our convention for which direction an edge points in.
  }

  // Stuff field pnext for each edge ie.
  for (auto ie=1; ie<=ne; ++ie) {
    auto aMin = 2.0 * M_PI;
    auto jeMin = -1;

    // Find the edge je incident to ihead(ie),
    // excluding the edge ie itself,
    // with smallest angle greater than ie's angle.
    // ("The Lowest Highest Point"! -- Moxy Fruvous)
    for (auto je=1; je<=ne; ++je) {
      if (je == ie)
	continue;
      if (ihead(ie) == ihead(je)) {
	const double a = Normalize(angle(je) - angle(ie));
	if (a < aMin) {
	  aMin = a;
	  jeMin = je;
	}
      }
      if (ihead(ie) == itail(je)) {
	const double a = Normalize(negangle(je) - angle(ie));
	if (a < aMin) {
	  aMin = a;
	  jeMin = -je;
	}
      }
    }
    rgedge[ie].pnext = jeMin;
    //printf("edge %d's pnext is %d\n", ie, jeMin);
  }

  // Stuff field mnext for each edge ie.
  for (auto ie=1; ie<=ne; ++ie) {
    auto aMin = 2.0 * M_PI;
    auto jeMin = -1;

    // Accumulate the edge je incident to itail(ie),
    // not ie itself,
    // with smallest angle greater than ie's angle.
    for (auto je=1; je<=ne; ++je) {
      if (je == ie)
	continue;
      if (itail(ie) == itail(je)) {
	const double a = Normalize(angle(je) - angle(ie));
	if (a < aMin) {
	  aMin = a;
	  jeMin = -je;
	}
      }
      if (itail(ie) == ihead(je)) {
	const double a = Normalize(negangle(je) - angle(ie));
	if (a < aMin) {
	  aMin = a;
	  jeMin = je;
	}
      }
    }
    rgedge[ie].mnext = jeMin;
  }

  // Stuff field ies for each vertex iv.
  for (auto iv=1; iv<=nv; ++iv) {
    // Find an edge incident to iv.
    for (auto ie=1; ie<=ne; ++ie) {
      if (ihead(ie) == iv) {
	rgvertex[iv].ies = ie;
	// iv'th vertex has an incident edge ie
	// (whose vertices are rgedge[ie].ihead and rgedge[ie].itail).
	break;
	}
      if (itail(ie) == iv) {
	rgvertex[iv].ies = -ie;
	// iv'th vertex has an incident edge ie
	// (whose vertices are rgedge[ie].itail and rgedge[ie].ihead).
	break;
      }
    }
    if (rgvertex[iv].ies == 0) {
      printf("Edahiro_Init got invalid input: %d'th vertex has no incident edge.\n", iv);
      return false;
    }
    //printf("vertex %d's ies is edge %d\n", iv, ies(iv));
  }

  return Preprocess();
}

// Returns -1 if point lies outside the hull.
int Edahiro_RegionFromPoint(double x, double y)
{
  return DoQuery(x, y) - 1;
}
