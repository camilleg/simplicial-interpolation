// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

// Simplicial interpolation.

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "si.h"
#include "bary.h"
#include "gacli.h"
#include "sammon.h"
#include "edahiro.h"

#include <GL/gl.h>
#include <GL/glu.h>
#ifdef __APPLE__
// Mac OS X 10.3.9 patch by Hans-Christoph Steiner, 2006 Jan 31
#include <glut.h>
#else
#include <GL/glut.h>
#endif

int cPoint = -1;
const char* szFileQi = "/tmp/.x";
const char* szFilePi = "/tmp/.y";
const char* szFileSi = "/tmp/.simplicial_interpolation.txt";

// Read d-vertices from szFile into an array.
// Caller is responsible for freeing memory.
vertex* readVertices(const char* szFile, int& count)
{
  // Count how many lines.
  count = 0;
  FILE* pf = fopen(szFile, "r");
  if (!pf)
    return NULL;
  char buf[1000];
  while (fgets(buf, sizeof(buf)-1, pf))
    ++count;
  rewind(pf);
  vertex* v = new vertex[count];
  for (int i=0; i<count; ++i)
    {
    if (fgets(buf, sizeof(buf)-1, pf) == NULL)
      fprintf(stderr, "Failed to read a line from '%s'.\n", szFile);
    char* pch = buf;
    char* pchNext;
    for (int j=0; j<d; ++j)
      {
      v[i][j] = strtod(pch, &pchNext);
      if (pchNext == pch)
        printf("Unexpected string \"%s\" while reading %dth float from %s\n",
	  pch, j, szFile);
      pch = pchNext;
      }
    }
  fclose(pf);
  return v;
}

// Read e-vertices from szFile into an array.
// Caller is responsible for freeing memory.
e_vertex* readEVertices(const char* szFile, int& count)
{
  // Count how many lines.
  count = 0;
  FILE* pf = fopen(szFile, "r");
  if (!pf)
    return NULL;
  char buf[1000];
  while (fgets(buf, sizeof(buf)-1, pf))
    ++count;
  rewind(pf);
  e_vertex* v = new e_vertex[count];
  for (int i=0; i<count; ++i)
    {
    if (fgets(buf, sizeof(buf)-1, pf) == NULL)
      fprintf(stderr, "Failed to read a line from e-vertices file '%s'.\n", szFile);
    char* pch = buf;
    char* pchNext;
    for (int j=0; j<e; ++j)
      {
      v[i][j] = strtod(pch, &pchNext);
      if (pchNext == pch)
        printf("Unexpected string \"%s\" while reading %dth float from %s\n",
	  pch, j, szFile);
      pch = pchNext;
      }
    }
  fclose(pf);
  return v;
}

// Read d-simplices from szFile into an array.
// count stores how many true simplices are read (no -1 vertex).
// countAll stores that, plus how many ray-simplices are read.
// True simplices precede ray-simplices in the returned array.
// Caller is responsible for freeing memory.
simplex* readSimplices(const char* szFile, int& count, int& countAll)
{
  // Count how many lines.
  count = 0;
  FILE* pf = fopen(szFile, "r");
  if (!pf)
    return NULL;
  char buf[1000];
  // First pass:  count lines.
  while (fgets(buf, sizeof(buf)-1, pf))
    ++countAll;
  simplex* sRet = new simplex[countAll];

  // Second pass:  collect true simplices.
  rewind(pf);
  int i1 = 0;
  for (int i=0; i<countAll; ++i)
    {
    simplex& s = sRet[i1];
    if (fgets(buf, sizeof(buf)-1, pf) == NULL)
      fprintf(stderr, "Failed to read a line from simplices file '%s'.\n", szFile);
    char* pch = buf;
    char* pchNext;
    bool fRay = false;
    for (int j=0; j<d+1; ++j)
      {
      s[j] = int(strtol(pch, &pchNext, 10));
      if (pchNext == pch)
        printf("Unexpected string \"%s\" while reading %dth int from %s\n",
	  pch, j, szFile);
      pch = pchNext;
      if (s[j] < 0)
        fRay = true;
      }
    if (!fRay)
      if (++i1 >= countAll)
        break;
    }
  count = i1;

  // Third pass:  collect ray-simplices.
  rewind(pf);
  for (int i=0; i<countAll; ++i)
    {
    simplex& s = sRet[i1];
    if (fgets(buf, sizeof(buf)-1, pf) == NULL)
      fprintf(stderr, "Failed to read a ray-simplex line from file '%s'.\n", szFile);
    char* pch = buf;
    char* pchNext;
    int jRay = -1;
    for (int j=0; j<d+1; ++j)
      {
      s[j] = int(strtol(pch, &pchNext, 10));
      pch = pchNext;
      if (s[j] < 0)
        jRay = j;
      }
    if (jRay >= 0)
      {
      // Disambiguate:  move -1 to the last position.
      s[jRay] = s[d];
      s[d] = -1;

      // Don't overflow sRet[].
      if (++i1 >= countAll)
        break;
      }
    }

  fclose(pf);
  if (i1 != countAll)
    printf("internal error in readSimplices.\n");
  return sRet;
}

vertex qC; // constructed common point of the ray-simplices
vertex* qi = NULL;
e_vertex* pi = NULL;
simplex* si = NULL;
simplexHint* hi = NULL;
int csi = 0;
int csiAll = 0;

bool dump_qi()
{
  // Give q_i values to hull, to prepare for Delaunay triangulation.
  FILE* pf = fopen(szFileQi, "w");
  if (!pf)
    return false;
  for (int i=0; i<cPoint; ++i)
    {
    for (int j=0; j<d; ++j)
      fprintf(pf, "%d ", int(qi[i][j]));
    fprintf(pf, "\n");
    }
  fclose(pf);
  return true;
}

bool init()
{
  if (e < d)
    {
    printf("error: e (%d) cannot be smaller than d (%d).\n", e, d);
    return false;
    }
  char buf[1000];

  // Get output sites p_i.

  const int cPointBecomesThis = d==2 ? 30 : d+10;
  sprintf(buf, "./rsites %d %d > %s", cPointBecomesThis, e, szFilePi);
  if (system(buf) < 0)
    fprintf(stderr, "system command failed\n");
  pi = readEVertices(szFilePi, cPoint);
  // szFilePi is a sequence of lines.
  // Each line is the coords of an e-vertex.
  // Vertex indices refer to this sequence (which starts at zero).

  // Get input sites q_i.

  enum { q_i_Manual, q_i_SammonsMapping, q_i_GeneticAlgorithm };
  const int q_i_kind = q_i_GeneticAlgorithm;
  switch (q_i_kind)
    {
  case q_i_Manual:
    sprintf(buf, "./rsites %d %d > %s", cPoint, d, szFileQi);
    if (system(buf) < 0)
      fprintf(stderr, "system command failed\n");
    qi = readVertices(szFileQi, cPoint);
    break;
  case q_i_SammonsMapping:
    qi = computeSammon(pi, cPoint, 1e6 /*scaled to match rsites's output*/);
    if (!dump_qi())
      return false;
    break;
  case q_i_GeneticAlgorithm:
    double* tmp = new double[e*cPoint];
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<e; ++j)
      tmp[i*e + j] = pi[i][j];
    const Member* foo = GADistanceMatrix(cPoint, e, d, tmp);
    delete [] tmp;
    qi = new vertex[cPoint];
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<d; ++j)
      qi[i][j] = double(foo->rgl[i*d + j]) / sHuge * 1e6;
      // *1e6 scales to match range of rsites's output.
    if (!dump_qi())
      return false;
    break;
    }
  if (!qi)
    return false;

  // Compute Delaunay triangulation.
  sprintf(buf, "cat %s | ./hull -d | sed '2,$!d' > %s", szFileQi, szFileSi);
  if (system(buf) < 0)
    fprintf(stderr, "system command failed\n");
  si = readSimplices(szFileSi, csi, csiAll);

  // printf("read %d true simplices, %d ray-simplices.\n", csi, csiAll-csi);

  // Output is a sequence of lines.  Each line lists the vertex indices of a d-simplex.
  // An index of -1 indicates the point at infinity;  we'll use that for
  // building ray-simplices.

  if (d == 2 && !Edahiro_Init(cPoint, qi, csi, si))
    return false;

  // Precompute some things to speed up eval().

  // Compute qC, centroid of bounding rectangle of all simplices.
  int i;
  for (i=0; i<d; ++i)
    {
    float zMin = DBL_MAX;
    float zMax = -DBL_MAX;
    for (int j=0; j<cPoint; ++j)
      {
      const double z = qi[j][i];
      if (z < zMin)
        zMin = z;
      if (z > zMax)
        zMax = z;
      }
    qC[i] = (zMin + zMax) * .5;
    }

  // Compute centroids of true simplices, and of the true simplices which are the core of the ray-simplices.
  vertex* viCentroid = new vertex[csiAll];
  for (i=0; i<csiAll; ++i)
    {
    vertex& v = viCentroid[i];
    const simplex& s = si[i];
    // Accumulate into v the centroid of s.
    for (int j=0; j<d; ++j)
      {
      // Accumulate into the j'th coord of v, the j'th coord of the centroid of s.
      // By adding the j'th coords of each vertex of s.
      double& t = v[j];
      t = 0.;
      for (int j1=0; j1<d+1; ++j1)
	{
	const int iVertex = s[j1];
	const vertex& v1 = iVertex<0 ? qC : qi[iVertex];
	t += v1[j];
	}
      t /= d+1;
      }
    }

  hi = new simplexHint[csiAll];
  for (i=0; i<csiAll; ++i)
    {
    if (!precomputeBary(si[i], hi[i], viCentroid[i], qi, &qC, i>=csi))
      return false;
    }

  delete [] viCentroid;
  return true;
}

void terminate()
{
  delete [] pi;
  delete [] qi;
  delete [] si;
}

// Return which member of the array of simplices si[] contains q.
// Fill w[0 to d+1] with q's barycentric coordinates w.r.t. that simplex.

int searchRaySimplices(const vertex& q, double* w)
{
  for (int i=csi; i<csiAll; ++i)
    if (computeBary(hi[i], q, w, true))
      return i;
  printf("internal error in searchRaySimplices\n");
  (void)computeBary(hi[0], q, w);
  return 0;
}

int searchBruteForce(const vertex& q, double* w)
{
  // Compute q's bary-coords w.r.t. each simplex s.
  // If one has coordinates all nonnegative, return that one.

  for (int i=0; i<csi; ++i)
    if (computeBary(hi[i], q, w))
      return i;
  return searchRaySimplices(q, w);
}

int searchEdahiro(const vertex& q, double* w)
{
  // Edahiro's algorithm handles only the case d==2.

  const int i = Edahiro_RegionFromPoint(q[0], q[1]);
  if (i >= 0)
    {
    if (i >= csi)
      {
      printf("internal error: edahiro returned out-of-range value\n");
      return 0;
      }
    if (!computeBary(hi[i], q, w))
      {
      printf("internal error: edahiro returned wrong simplex.\n");
      return i;
      }
    return i;
    }
  return searchRaySimplices(q, w);
}

#define DATAVIZ
#ifdef DATAVIZ
int iFound = 0;
double wFound[d+1] = {0};
vertex rFound;
#endif

// Map d-vertex q to e-vertex p.
void eval(const vertex& q, e_vertex& p)
{
  // Find which simplex s contains q.

  double w[d+1]; // q's coordinates w_j with respect to s.
  const int iS = d==2 ? searchEdahiro(q, w) : searchBruteForce(q, w);
  const simplex& s = si[iS];
  int i, j;

  //printf("Barycoords ");
  //for (i=0; i<d+1; ++i)
  //  printf("%g  ", w[i]);
  //printf("\n");

#ifdef TESTING
  // Verify that q == the point whose barycoords are w[] wrt s.
  vertex r;
  for (j=0; j<d; ++j)
    {
    r[j] = 0.;
    for (i=0; i<d+1; ++i)
      {
      const vertex& v = s[i] < 0 ? qC : qi[s[i]];
      r[j] += w[i] * v[j];
      }
    }
  // Now r is the reconstructed point.  How far is it from q?
  // (i.e., how accurate were the barycentric coordinates w[]?
  double dist = 0.;
  for (j=0; j<d; ++j)
    dist += sq(r[j] - q[j]);
  dist = sqrt(dist);
  if (dist > 1e-8)
    printf("warning: reconstruction error = %g\n\n\n", dist);
#ifdef DATAVIZ
  rFound = r;
#endif
#endif

#ifdef DATAVIZ
  iFound = iS;
  for (j=0; j<d+1; ++j)
    wFound[j] = w[j];
#endif

  // Compute a weighted sum, weights w[] and vertex pi[s[]].

  for (j=0; j<e; ++j)
    {
    double _ = 0;
    for (i=0; i<d+1; ++i)
      {
      _ += w[i] * pi[s[i]][j];
      }
    p[j] = _;
    }
}

void evalAutomatic()
{
  // Exercise the evaluation function.
  char buf[1000];
  const char* szFileTest = "/tmp/.test";
  const int ctest = 20; // 100000 for timing tests
  sprintf(buf, "./rsites %d %d > %s", ctest, d, szFileTest);
  if (system(buf) < 0)
    fprintf(stderr, "system command failed\n");

  // Compiled -O2 or -O3 for a 1GHz Pentium III,
  // this does one test in 310 usec or 3200 tests per second for e=42, d=7.
  // For d=3, 100 usec/test or 10000 tests/sec.
  // For d=2, 66 usec/test or 15000 tests/sec.

  int cqtest;
  vertex* qtest = readVertices(szFileTest, cqtest);
  for (int i=0; i<cqtest; ++i)
    {
    printf("eval %d/%d\n", i, cqtest);
    e_vertex p;
    qtest[i].dump("query: ");
    eval(qtest[i], p);
    p.dump("result: ");
    }
}

void drawCircle(double radius)
{
  // optimize: precompute cos&sin;  scale instead of *radius.
  glBegin(GL_POLYGON);
  for (float a=0; a<=2*M_PI; a += .1)
    glVertex2f(radius * cos(a), radius * sin(a));
  glEnd();
}

vertex vQ;

void display()
{
  glClear(GL_COLOR_BUFFER_BIT);
  int i, j;

  const bool fInside = iFound < csi;
  const simplex& s = si[iFound];

  // Bounding simplex (i.e., triangle).
  glColor3f(.05,.05,.4);
  glBegin(GL_TRIANGLES);
  if (fInside)
    {
    for (j=0; j<=d; ++j)
      {
      const vertex& v = qi[s[j]];
      glVertex2f(v[0], v[1]);
      }
    }
  else
    {
    glVertex2f(qC[0], qC[1]);
    vertex v = qi[s[0]];
    v[0] += 100. * (v[0] - qC[0]);
    v[1] += 100. * (v[1] - qC[1]);
    glVertex2f(v[0], v[1]);
    v = qi[s[1]];
    v[0] += 100. * (v[0] - qC[0]);
    v[1] += 100. * (v[1] - qC[1]);
    glVertex2f(v[0], v[1]);
    }
  glEnd();

  // Disc around each vertex of containing simplex.
  // Disc's radius indicates weight (barycentric coordinate) of corresponding vertex.
  glColor3f(.5,.5,.5);
  for (i=0; i<=d; ++i)
    {
    glPushMatrix();
    const int iVertex = s[i];
    const vertex& v = iVertex<0 ? qC : qi[iVertex];
    glTranslatef(v[0], v[1], 0);
    drawCircle(sqrt(fabs(wFound[i])) * 5e4);
    glPopMatrix();
    }

  // Triangulation mesh.
  if (fInside)
    glColor3f(1,0,.4);
  else
    glColor3f(.5,0,.2);
  for (i=0; i<csi; ++i)
    {
    const simplex& s = si[i];
    glBegin(GL_LINE_LOOP);
    for (j=0; j<=d; ++j)
      {
      const vertex& v = qi[s[j]];
      glVertex2f(v[0], v[1]);
      }
    glEnd();
    }

  // Ray-simplices.
  if (fInside)
    glColor3f(0,.5,0);
  else
    glColor3f(0,1,0);
  for (i=csi; i<csiAll; ++i)
    {
    const simplex& s = si[i];
    glBegin(GL_LINE_LOOP);
      glVertex2f(qC[0], qC[1]);
      const vertex& v = qi[s[0]]; glVertex2f(v[0], v[1]);
      const vertex& w = qi[s[1]]; glVertex2f(w[0], w[1]);
    glEnd();
    }

  glColor3f(1,1,1);
  glRasterPos2f(qC[0], qC[1]);
  glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'C');

#ifdef TESTING
  // Reconstructed point.
  glColor3f(1,.2,.2);
  glRasterPos2f(rFound[0], rFound[1]);
  glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'R');
#endif

  // Query point.
  glColor3f(1,1,1);
  glRasterPos2f(vQ[0], vQ[1]);
  glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'q');

  // Vertices.
  glColor3f(1,1,0);
  for (i=0; i<cPoint; ++i)
    {
    glRasterPos2f(qi[i][0], qi[i][1]);
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '0' + i);
    }

  glutSwapBuffers();
  usleep(2000); // CPU throttle
}

int xSize = 700;
int ySize = 700;

const int margin = 200000;

inline void XYFromMouse(double& x, double& y, int xM, int yM)
{
  x = double(xM) / double(xSize)        * (1e6 + 2*margin) - margin;
  y = (1. - double(yM) / double(ySize)) * (1e6 + 2*margin) - margin;
}

void mouse(int button, int state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
    XYFromMouse(vQ[0], vQ[1], x, y);
    e_vertex p;
    eval(vQ, p);
    }
}

void drag(int x, int y)
{
  mouse(GLUT_LEFT_BUTTON, GLUT_DOWN, x, y);
}

void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
  switch (key)
    {
  case 27: /* escape key */
    terminate();
    exit(0);
    }
}

void reshape(int w, int h)
{
  glViewport(0, 0, w, h);
  xSize = w;
  ySize = h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-margin, margin+1e6, -margin, margin+1e6);
}

void evalInteractive(int argc, char** argv)
{
  // Evaluate points interactively.
  // (For verifying the barycentric coords code.)

  if (d != 2)
    {
    printf("oops, evalInteractive requires that d equals 2.\n");
    return;
    }

  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(xSize,ySize);
  glutCreateWindow("Simplicial Interpolator");
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(drag);
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutIdleFunc(display);
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();
  glutMainLoop(); // never returns
}

int main(int argc, char** argv)
{
  if (!init())
    return -1;
  // evalAutomatic();
  evalInteractive(argc, argv);
  terminate();
  return 0;
}
