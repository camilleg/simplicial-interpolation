// Approximate a solution with a genetic algorithm.
// Then relax to the local maximum that the GA was approaching.

#include <algorithm>
#include <cmath>
#include <cstdlib> // random()
#include <iostream>
#include <limits>

#undef VERBOSE
#ifdef VERBOSE
#include <cstdio>
#endif

#include "ga.h"
#include "gacli.h"
#include "util.h"

// All of this file's vars and functions could be wrapped into a class, of course.
static int cpt = -1;
static int cdimSrc = -1;
static int cdimDst = -1;

//;; iloop jloop rg[triangularNumber(cpt,i,j)]  -->  iloop rg[i]

template <class T> T sq(const T& _) { return _*_; }

/*
Each linear dimension of the space may have a vastly different
range (in the millions or in the millionths) from others,
so normalize it somehow.  As a first guess, assume all dimensions
are equally important:  scale the max difference along each dimension
to be the same (if nonzero).
Later: manual scaling overrides (multiply rgzScale[idim] by that).
Later: log-scale distance option for some dimensions. (use log(x) internally
       in the linear geometry; use x for i/o.
*/

void InitDistanceMatrixZ(int cpt, int cdimSrc, double* rgzDist, const double* rgzPt)
{
  // compute scaling factors for each dimension
  double rgzScale[cdimSrc];
  for (auto idim=0; idim<cdimSrc; ++idim) {
    auto zMin = std::numeric_limits<double>::max();
    auto zMax = -zMin;
    for (auto i=0; i<cpt; ++i) {
      zMin = std::min(zMin, rgzPt[i * cdimSrc + idim]);
      zMax = std::max(zMax, rgzPt[i * cdimSrc + idim]);
    }
    rgzScale[idim] = zMax==zMin ? 1.0 : 1.0/(zMax-zMin);
  }

  auto zDistMax = 0.0;
  for (auto i=0; i<cpt-1; ++i)
  for (auto j=i+1; j<cpt; ++j) {
    auto zSum = 0.0;
    for (auto idim = 0; idim < cdimSrc; ++idim) {
      zSum += sq((rgzPt[i * cdimSrc + idim] -
	          rgzPt[j * cdimSrc + idim]) * rgzScale[idim]);
    }
    zDistMax = std::max(zDistMax, rgzDist[TriIJ(i,j,cpt)] = sqrt(zSum));
  }
  if (zDistMax == 0.0) {
    std::cerr << "InitDistanceMatrixZ() got inputs that made only zero distances.  Crash imminent.\n";
    // rgzPt[] might have been all zeros.  Just DBZ and get it over with.
  }
  // Normalize distances wrt longest distance.
  for (auto i=0; i<cpt-1; ++i)
  for (auto j=i+1; j<cpt; ++j)
    rgzDist[TriIJ(i,j,cpt)] /= zDistMax;
}

void InitDistanceMatrixL(int cpt, int cdimDst, double* rgzDist, short* rgzPt)
{
  auto zDistMax = 0.0;
  for (auto i=0; i<cpt-1; ++i)
  for (auto j=i+1; j<cpt; ++j) {
    short* psi = rgzPt + i*cdimDst;
    short* psj = rgzPt + j*cdimDst;
    auto zSum = 0.0;
    for (auto idim = cdimDst; idim > 0; --idim)
      zSum += sq(*psi++ - *psj++);
    zDistMax = std::max(zDistMax, rgzDist[TriIJ(i,j,cpt)] = sqrt(zSum));
  }
  // normalize distances wrt longest distance.
//;; this can be a k-loop, and above in InitDistanceMatrixZ
  for (auto i=0; i<cpt-1; ++i)
  for (auto j=i+1; j<cpt; ++j)
    rgzDist[TriIJ(i,j,cpt)] /= zDistMax;
}

// Compute RMS error between 2 distance vectors (possibly triangular matrices).
inline double DDistanceMatrix(double* rgzDist0, double* rgzDist1, int cpt)
{
#ifdef VERBOSE
  for (auto k=0; k<cpt; ++k) {
    if (std::isnan(rgzDist0[k])) {
      std::cerr << "DDistanceMatrix got corrupt input.";
      return std::numeric_limits<double>::signaling_NaN();
    }
  }
#endif
  auto z = 0.0;
  for (auto k=0; k<cpt; ++k)
    z += sq(rgzDist0[k] - rgzDist1[k]);

//return sqrt(z / (double)cpt);
//sqrt and 1/cpt are monotonic and constant over the GA's execution,
//so don't bother with them.
  return z;
}

constexpr auto zTweakBuffer = 0.8; // chicken factor

// scale the member up, same scale along all dimensions.
void Tweak(void* pv)
{
  short sMin = sHuge;
  short sMax = -sHuge;
  short* ps = ((Member*)pv)->rgl;
  for (auto i=0; i<cpt*cdimDst; ++i) {
    sMin = std::min(sMin, ps[i]);
    sMax = std::max(sMax, ps[i]);
  }

  // Linearly map [sMin, sMax] to [0, sHuge].
  //;; Skip this if m is almost 1?  Collect stats on value of m.
  const auto m = zTweakBuffer * sHuge / (sMax - sMin);
  for (auto i=0; i<cpt*cdimDst; ++i)
    ps[i] = (ps[i] - sMin) * m;
}

void GenerateRandom(void* pv)
{
  short* ps = ((Member*)pv)->rgl;
  for (auto i=0; i<cpt*cdimDst; ++i)
    ps[i] = short(random() & sHuge);
  Tweak(pv);
}

void MutateRandom(void* pv, long cIter)
{
  static long lMask = random();
  {
    // Use the same bit pattern for a while,
    // so the same points get tweaked in one generation.
    static int _= -1;
    if ((++_ %= 50) == 0) lMask = random();
  }

  auto ibit = 0;
  const int cbit = 2; // Size of bitfield to test.  Bigger means fewer mutations.
  const short denom = 1.0 + 3.0 * sqrt(10.0 * cIter);
  short* ps = ((Member*)pv)->rgl;
  for (auto i=0; i<cpt*cdimDst; ++i) {
    // Are all cbit bits of the mask, starting at the ibit'th, zero?
    if ((lMask & (long((1 << cbit) - 1) << ibit)) == 0)
      ps[i] += ((random() & sHuge) - sHuge/2) / denom;
    // Try next cbit bits next time.
    if ((ibit += cbit) >= 8 * int(sizeof(long)) - cbit)
      ibit = 0;
  }
  Tweak(pv);
}

static double* rgzDistSrc = nullptr;
static double* rgzDistDst = nullptr;

double ComputeSuitability(void* pv)
{
  InitDistanceMatrixL(cpt, cdimDst, rgzDistDst, ((Member*)pv)->rgl);
  return -DDistanceMatrix(rgzDistSrc, rgzDistDst, triangularNumber(cpt-1));
  // Larger distance means less suitable.
}

// Caller free()s the return value.
Member* GADistanceMatrix(int cptArg, int cdimSrcArg, int cdimDstArg, double* rgzPt)
{
  cpt = cptArg;
  cdimSrc = cdimSrcArg;
  cdimDst = cdimDstArg;

  // If these two could be passed to ComputeSuitability(), they could be local vars.
  rgzDistSrc = new double[triangularNumber(cpt-1)];
  rgzDistDst = new double[triangularNumber(cpt-1)];

  InitDistanceMatrixZ(cpt, cdimSrc, rgzDistSrc, rgzPt);

  const auto pmemberBest = (Member*)GA(
    sizeof(short) * cpt * cdimDst,
    GenerateRandom,
    MutateRandom,
    Tweak,
    ComputeSuitability,
    0.0,
    50,
    0.6 /* timeout, in seconds */
  );

  for (auto i=0; i<cpt*cdimDst; ++i)
    pmemberBest->rgl[i] /= zTweakBuffer;

#ifdef VERBOSE
  // Print the matrix.
  InitDistanceMatrixL(cpt, cdimDst, rgzDistDst, pmemberBest->rgl);
  for (auto i=0; i<cpt-1; ++i)
    {
    for (auto j=i+1; j<cpt; ++j)
      printf("%.02f  ", rgzDistDst[TriIJ(i,j,cpt)]);
    putchar('\n');
    }
#endif
  delete [] rgzDistDst;
  delete [] rgzDistSrc;
  return pmemberBest;
}
