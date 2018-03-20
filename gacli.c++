/************************************************************************

This library implements simplicial interpolation as described in
"Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
Cambridge University Press.

Copyright 2002 Camille Goudeseune.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation;  either version 2.1
of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Contact:
  Camille Goudeseune
  Integrated Systems Laboratory
  Beckman Institute for Advanced Science and Technology
  405 N Mathews
  Urbana IL 61801 USA
  camilleg@isl.uiuc.edu

************************************************************************/

// Use GA to get an approximate solution,
// then use relaxation/iteration to converge to the local maximum which
// the GA was approaching.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "ga.h"
#include "gacli.h"
#include "util.h"

// Make these class variables if you want reentrancy.
static int cpt = -1;
static int cdimSrc = -1;
static int cdimDst = -1;

//;; iloop jloop rg[triangularNumber(cpt,i,j)]  -->  iloop rg[i]

double* rgzDistSrc;
double* rgzDistDst;

inline double sq(double _) { return _*_; }

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

void InitDistanceMatrixZ(int cpt, int cdimSrc, double* rgzDist, double* rgzPt)
{
  int i;
  int j;
  double rgzScale[100/*;;max dim*/];

  {
  // compute scaling factors for each dimension
    for (int idim=0; idim<cdimSrc; idim++)
      {
      double zMin =  DBL_MAX;
      double zMax = -DBL_MAX;
      for (i=0; i<cpt; i++)
	{
	if (rgzPt[i * cdimSrc + idim] < zMin)
	  zMin = rgzPt[i * cdimSrc + idim];
	if (rgzPt[i * cdimSrc + idim] > zMax)
	  zMax = rgzPt[i * cdimSrc + idim];
	}
      rgzScale[idim] = (zMax==zMin) ? 1.f : 1.f / (zMax - zMin);
      }
  }

  double zDistMax = 0.;
  for (i=0; i<cpt-1; i++)
  for (j=i+1; j<cpt; j++)
    {
    double zSum = 0.;
    for (int idim = 0; idim < cdimSrc; idim++)
      {
      zSum += sq((rgzPt[i * cdimSrc + idim] -
	 rgzPt[j * cdimSrc + idim]) * rgzScale[idim]);
      }
    double tmp = rgzDist[TriIJ(i,j,cpt)] = sqrt(zSum);
    if (tmp > zDistMax)
      zDistMax = tmp;
    }
  // normalize distances wrt longest distance.
  for (i=0; i<cpt-1; i++)
  for (j=i+1; j<cpt; j++)
    rgzDist[TriIJ(i,j,cpt)] /= zDistMax;
}

void InitDistanceMatrixL(int cpt, int cdimDst, double* rgzDist, short* rgzPt)
{
  int i;
  int j;

  double zDistMax = 0.;
  for (i=0; i<cpt-1; i++)
  for (j=i+1; j<cpt; j++)
    {
    short* psi = rgzPt + i*cdimDst;
    short* psj = rgzPt + j*cdimDst;
    double zSum = 0.;
    for (int idim = cdimDst; idim > 0; idim--)
      zSum += sq((double)(*psi++ - *psj++));
    double tmp = rgzDist[TriIJ(i,j,cpt)] = sqrt(zSum);
    if (tmp > zDistMax)
      zDistMax = tmp;
    }
  // normalize distances wrt longest distance.
//;; this can be a k-loop, and above in InitDistanceMatrixZ
  for (i=0; i<cpt-1; i++)
    {
    for (j=i+1; j<cpt; j++)
      rgzDist[TriIJ(i,j,cpt)] /= zDistMax;
    }
}

// Compute RMS error between 2 distance vectors (possibly triangular matrices).
inline double DDistanceMatrix(double* rgzDist0, double* rgzDist1, int cpt)
{
  double z = 0.;
  for (int k=0; k<cpt; ++k)
    z += sq(rgzDist0[k] - rgzDist1[k]);

//return sqrt(z / (double)cpt);
//sqrt and 1/cpt are monotonic and constant over the GA's execution,
//so don't bother with them.
  return z;
}

const double zTweakBuffer = .8f; // chicken factor

// scale the member up, same scale along all dimensions.
void Tweak(void* pv)
{
  short sMin = sHuge;
  short sMax = -sHuge;

  int i;
  short* ps = ((Member*)pv)->rgl;
  for (i=0; i<cpt*cdimDst; ++i)
    {
    if (ps[i] < sMin)
      sMin = ps[i];
    if (ps[i] > sMax)
      sMax = ps[i];
    }

  // now linearly map [sMin, sMax] to [0, sHuge]
  double m = zTweakBuffer * (double)sHuge / (double)(sMax - sMin);

  //;; skip this if m is almost 1?  Collect stats on value of m.
  for (i=0; i<cpt*cdimDst; ++i)
    ps[i] = (short)((ps[i] - sMin) * m);
}

void GenerateRandom(void* pv)
{
  short* ps = ((Member*)pv)->rgl;
  for (int i=0; i<cpt*cdimDst; ++i)
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

  int ibit = 0;
  const int cbit = 2;	// size of bitfield to test.
			// big == fewer mutations.
  const short denom = short(1 + 3. * sqrt(10.*cIter));
  short* ps = ((Member*)pv)->rgl;
  for (int i=0; i<cpt*cdimDst; ++i)
    {
    // Are all cbit bits of the mask, starting at the ibit'th, zero?
    if ((lMask & (((1 << cbit) - 1) << ibit)) == 0)
      ps[i] += short(((random() & sHuge) - sHuge/2) / denom);
    // Try next cbit bits next time.
    if ((ibit+=cbit) >= int(sizeof(long))*8-cbit)
      ibit = 0;
    }
  Tweak(pv);
}

double ComputeSuitability(void* pv)
{
  InitDistanceMatrixL(cpt, cdimDst, rgzDistDst, ((Member*)pv)->rgl);
  return -DDistanceMatrix(rgzDistSrc, rgzDistDst, triangularNumber(cpt-1));
  // Larger distance == less suitable.
}

Member* GADistanceMatrix(int cptArg, int cdimSrcArg, int cdimDstArg, double* rgzSrc)
{
  cpt = cptArg;
  cdimSrc = cdimSrcArg;
  cdimDst = cdimDstArg;

  rgzDistSrc = (double *)calloc(triangularNumber(cpt-1), sizeof(double));
  rgzDistDst = (double *)calloc(triangularNumber(cpt-1), sizeof(double));

  InitDistanceMatrixZ(cpt, cdimSrc, rgzDistSrc, rgzSrc);

  static Member* pmemberBest;

  pmemberBest = (Member*)GA(
    sizeof(short) * cpt * cdimDst,
    GenerateRandom,
    MutateRandom,
    Tweak,
    ComputeSuitability,
    0.,
    50,
    0.6 /* timeout, in seconds */
    );

  {
  short* ps = pmemberBest->rgl;
  for (int i=0; i<cpt*cdimDst; i++)
    ps[i] = (short)((double)(ps[i]) / zTweakBuffer);
  }

#ifdef VERBOSE
  {
  // just to print out matrix
  InitDistanceMatrixL(cpt, cdimDst, rgzDistDst, pmemberBest->rgl);
  for (int i=0; i<cpt-1; i++)
    {
    for (int j=i+1; j<cpt; j++)
      printf("%.3g  ", rgzDistDst[TriIJ(i,j,cpt)]);
    putchar('\n');
    }
  }
#endif
  free(rgzDistSrc);
  free(rgzDistDst);
  return pmemberBest;
}
