#include <cmath>
#include <cstdio>
#include <cstdlib> // drand48(), rand()
#include <cstring>
#include <limits>

#include "sammon.h"
#include "util.h"

void computeSammon(std::vector<vertex>& qi, const std::vector<vertex>& pi, const double scalar)
{
  // Modification of the published algorithm:  store the squares of
  // the distances, instead of the distances themselves.
  // Avoids computing a zillion square roots.

  const int d = qi.front().size();
  const int e = pi.front().size();

  // Compute target values for distance matrix.
  const int cpt = pi.size();
  double rgzTarget[triangularNumber(cpt - 1)];
  int i, j, k=0;
  for (i=0; i<cpt-1; ++i)
  for (j=i+1; j<cpt; ++j, ++k)
    {
    rgzTarget[k] = 0.;
    for (int idim=0; idim<e; ++idim)
      rgzTarget[k] += sq(pi[i][idim] - pi[j][idim]);
    }

  // Run Sammon's Mapping crun times, and keep the best run so far.
  double rgz[cpt*d];
  double rgzBest[cpt*d];
  const int crun = 1000;
  int runLastGood = -1;
  double errorMin = std::numeric_limits<double>::max();
  for (int run=0; run<crun; ++run)
    {
    if (run - runLastGood > crun/4)
      // Haven't improved in quite a while.
      break;

    // Build initial configuration of points.
    for (i=0; i<cpt*d; ++i)
      rgz[i] = drand48() * scalar;
    
    // Now rgz[i*d + 0] to rgz[i*d + d-1] are the coords of
    // the i'th point of the configuration.

    // Iterate.
    const int iterMax = 75 * cpt;
    for (int iter=0; iter<iterMax; iter++)
      {
      // Pick a pair of points uniformly randomly (fast and simple).
      do {
	i = (rand() >> 4) % cpt;
	j = (rand() >> 4) % cpt;
	}
      while (i == j);
      if (i>j) std::swap(i, j);
      // Now 0 <= i < j < cpt.

      // Find the distance between them, target and current.
      const double distTarget = rgzTarget[TriIJ(i, j, cpt)];

      double* aa = &rgz[i*d];
      double* bb = &rgz[j*d];
      const double* a = aa;
      const double* b = bb;
      double distCur = 0.;
      for (k=0; k<d; ++k)
	distCur += sq(b[k] - a[k]);

      // Adjust their relative positions.
      const double temperature = 1. - (double)iter/iterMax;     // From 1 to 0.
      const double gamma = .8 * temperature * temperature;      // Also down to 0.
      const double magnitude = gamma * (distTarget - distCur) / distCur;
      for (k=0; k<d; ++k)
	{
	const double c = magnitude * (a[k]-b[k]);
	aa[k] += c;
	bb[k] -= c;
	}
      
      // Now rgz's pairwise distances approximate rgzTarget's.
      // See how good the approximation really is.
      // (Later optimization:  update error only incrementally, since only two points changed.)

      double error = 0.;
      k = 0;
      for (i=0; i<cpt-1; ++i)
      for (j=i+1; j<cpt; ++j, ++k)
	{
	const double* a = &rgz[i*d];
	const double* b = &rgz[j*d];
	double sum = 0.;
	for (int _=0; _<d; ++_)
	  sum += sq(b[_] - a[_]);
	error += fabs(rgzTarget[k] - sum);
	}
      if (error < errorMin)
	{
	// Best approximation so far.  Keep this one.
	// (Later optimization:  reduce gamma and hang around here for a while,
	// lest we wander away from the local minimum near here.)
	errorMin = error;
	runLastGood = run;
	memcpy(rgzBest, rgz, cpt*d*sizeof(double));
	//if (run > 0)
	//  printf("\t%5d/%.3f: %.3f\n", run, float(iter)/iterMax, sqrt(error));
	}
      }
    }

  // printf("Sammon's mapping: best error is %.3f\n", sqrt(errorMin));
  for (i=0; i<cpt; ++i)
  for (k=0; k<d; ++k)
    qi[i][k] = rgzBest[i*d + k];
}
