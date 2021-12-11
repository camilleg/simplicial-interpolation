// Derived closely from vss's map/mymain.c++,
// which derives from hullmain.c++.
// Uses many global variables from the Hull code.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <random>

#include "hull.h"
#include "si.h"

static point site_blocks[MAXBLOCKS];
static int num_blocks = 0;
static long num_sites = 0L;
static int dim = -1;
FILE* DFILE = stderr;

int vcpt = 0;
HH* vpH = nullptr;
TT* vpT = nullptr;
int iH = 0;
int iT = 0;

long site_numm(site p) {
  if (p==hull_infinity)
    return -1;
  if (!p)
    return -2;
  long j;
  for (long i=0L; i<num_blocks; ++i)
    if ((j=p-site_blocks[i])>=0 && j < BLOCKSIZE*dim) 
      return j/dim + BLOCKSIZE*i;
  return -3;
}

site new_site(site p, long j) {
  assert(num_blocks+1 < MAXBLOCKS);
  if (0 == j%BLOCKSIZE) {
    assert(num_blocks < MAXBLOCKS);
    return site_blocks[num_blocks++] = (site)malloc(BLOCKSIZE*site_size);
  }
  return p+dim;
}

site read_next_site(long j) {
  if (j >= vcpt)
    return nullptr; // end of list
  hull_p = new_site(hull_p,j);
  for (auto i=0; i<dim; ++i) {
    extern std::vector<vertex> qi;
    hull_p[i] = floor(mult_up * qi[j][i] + 0.5);    // these are the input points.
    if (hull_p[i] < mins[i])
      mins[i] = hull_p[i];
    if (hull_p[i] > maxs[i])
      maxs[i] = hull_p[i];
  }
  return hull_p;
}

site get_site_offline(long i)
{
  return i>=num_sites ? nullptr : site_blocks[i/BLOCKSIZE] + (i%BLOCKSIZE)*dim;
}

std::vector<long> shufmat;
void make_shuffle() {
  shufmat.resize(num_sites + 1);
  std::iota(shufmat.begin(), shufmat.end(), 0);
  std::shuffle(shufmat.begin(), shufmat.end()-1, std::default_random_engine{std::random_device{}()});
}

static long s_num = 0L;
site get_next_site() {
  return get_site_offline(shufmat[s_num++]);
}

void make_output(
  simplex* root,
  simplex* (*visit_gen)(simplex*, visit_func*),
  visit_func visit,
  out_func* out_funcp,
  FILE *F)
{
  out_funcp(0,0,F,-1);
  visit(0, (void*)out_funcp);
  visit_gen(root, visit);
  out_funcp(0,0,F,1);
}

void resetEverything()
{
  extern long pnum;
  pnum = 0;
  num_blocks = 0;
  num_sites = 0;
  vpH = nullptr;
  vpT = nullptr;
  s_num = 0L;
  hull_p = nullptr;
  rdim = cdim = 0;
  check_overshoot_f = 0;
  set_ch_root(nullptr);
  mult_up = 1.0; // Already scaled by si.c++'s constexpr auto scale = 1e6.
}

bool delaunay_tri(std::vector<d_simplex>& si, std::vector<d_simplex>& siRay, int dimArg, int cPt) {
  resetEverything();
  if (cPt <= 1) {
    // If cPt==1, iH==iT==1 but vpH and vpT remain uninitialized.
    std::cerr << "Triangulation needs at least 2 points, not " << cPt << ".\n";
    goto LAbort;
  }
  if (dimArg > MAXDIM) {
    std::cerr << "dimension bound MAXDIM exceeded\n";
LAbort:
    si.clear();
    siRay.clear();
    return false;
  }
  dim = dimArg;
  point_size = site_size = dim * sizeof(Coord);

  // Read vcpt points from qi.
  vcpt = cPt;
  for (num_sites=0; read_next_site(num_sites); ++num_sites);

  make_shuffle();
  const auto root = build_convex_hull(get_next_site, site_numm, dim, 1);

  // Approximation.  CG_vlist_out stuffs vpT and vpH.  Buffer overflow.
  const auto maxTets = 100 * cPt * 4 + 3;

  // If vpT and vpH were on the stack, if maxTets > 200000 or so, the stack might overflow.
  vpT = new TT[maxTets]; // output: list of vertex indices, defining tetrahedra.
  vpH = new HH[maxTets]; // output: list of vertex indices, defining triangles on convex hull.
  iH = 0;
  iT = 0;

  // Stuff T[i], T[i]+1, T[i]+2, T[i]+3 with tetrahedra.
  // Stuff H[i], H[i]+1, H[i]+2         with triangles on the hull.
  make_output(root, visit_hull, facets_print, &CG_vlist_out, stdout);
  // CG_vlist_out wrote output to vpT,iT (simplices aka tets) and vpH,iH (ray-simplices aka hull, without the -1).

  extern std::vector<vertex> qi;
  const auto d = qi.front().size();
  resize(si,    d+1, iT); // aka ctet
  resize(siRay, d+1, iH); // aka ctri
  for (auto i=0; i<iT; ++i) {
    const auto& first = vpT[i].begin();
    std::copy(first, first + d+1, si[i].begin());
  }
  for (auto i=0; i<iH; ++i) {
    auto& s = siRay[i];
    // vpH[][d] is left uninitialized by CG_vlist_out().  It's implicitly -1.
    const auto& first = vpH[i].begin();
    std::copy(first, first + d, s.begin());
    s[d] = -1;
  }
  delete [] vpH;
  delete [] vpT;

  free_hull_storage();
  for (auto i=0; i<num_blocks; ++i)
    free(site_blocks[i]);
  return true;
}
