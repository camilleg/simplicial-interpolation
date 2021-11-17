// Simplicial interpolation.

#include <cmath>
#include <cstdio>
#include <cstdlib> // strtoXXX, erand48, exit
#include <limits>
#include <unistd.h>

#include "bary.h"
#include "edahiro.h"
#include "gacli.h"
#include "sammon.h"

#ifdef __APPLE__
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

int cPoint = -1;
const char* szFileQi = "/tmp/.x";
const char* szFileSi = "/tmp/.simplicial_interpolation.txt";

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

vertex qC{0}; // Constructed common point of the ray-simplices.
e_vertex pC{0}; // What qC maps to.
vertex* qi = nullptr;
e_vertex* pi = nullptr;
simplex* si = nullptr;
simplexHint* hi = nullptr;
auto csi = 0;
auto csiAll = 0;

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

// Scale the inputs to the hull algorithm, which uses exact integer arithmetic.
// Outside [1e2, 1e7], ch.c++ suffers degeneracies and overshoots.
constexpr auto scale = 1e6;

void randomSites_d(vertex* dst, int n) {
  unsigned short x[3]; // deliberately uninitialized for a mere demo
  for (auto j=0; j<n; ++j)
  for (auto i=0; i<d; ++i)
    dst[j][i] = scale * erand48(x);
}
void randomSites_e(e_vertex* dst, int n) {
  unsigned short x[3]; // deliberately uninitialized for a mere demo
  for (auto j=0; j<n; ++j)
  for (auto i=0; i<e; ++i)
    dst[j][i] = scale * erand48(x);
}

e_vertex eval(const vertex&);

bool init()
{
  if (e < d)
    {
    printf("error: e (%d) cannot be smaller than d (%d).\n", e, d);
    return false;
    }

  // Make output sites p_i.
  {
    const auto cPointBecomesThis = d==2 ? 30 : d+10;
    pi = new e_vertex[cPointBecomesThis];
    randomSites_e(pi, cPointBecomesThis);
    cPoint = cPointBecomesThis;
  }

  // Get input sites q_i.

  enum { q_i_Manual, q_i_SammonsMapping, q_i_GeneticAlgorithm };
  const int q_i_kind = q_i_GeneticAlgorithm;
  switch (q_i_kind)
    {
  case q_i_Manual:
    qi = new vertex[cPoint];
    randomSites_d(qi, cPoint);
    break;
  case q_i_SammonsMapping:
    qi = computeSammon(pi, cPoint, scale);
    if (!dump_qi())
      return false;
    break;
  case q_i_GeneticAlgorithm:
    double tmp[e*cPoint];
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<e; ++j)
      tmp[i*e + j] = pi[i][j];
    const Member* foo = GADistanceMatrix(cPoint, e, d, tmp);
    qi = new vertex[cPoint];
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<d; ++j)
      qi[i][j] = scale * double(foo->rgl[i*d + j]) / sHuge;
    if (!dump_qi())
      return false;
    break;
    }
  if (!qi)
    return false;

  // Compute Delaunay triangulation.
  char buf[1000];
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

  // Calculate qC, the centroid of the bounding box of all simplices.
  for (auto i=0; i<d; ++i)
    {
    auto zMin = std::numeric_limits<double>::max();
    auto zMax = -zMin;
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

  // Calculate the centroid of each true simplex,
  // and of each true simplex that is the core of a ray-simplex.
  vertex viCentroid[csiAll];
  for (auto i=0; i<csiAll; ++i)
    {
    auto& v = viCentroid[i];
    const simplex& s = si[i];
    // Accumulate into v the centroid of s.
    for (auto j=0; j<d; ++j)
      {
      // Set the j'th coord of v, the centroid of s,
      // to the mean of the j'th coords of s's vertices.
      v[j] = 0.0;
      for (auto j1=0; j1<d+1; ++j1)
	v[j] += (s[j1] < 0 ? qC : qi[s[j1]])[j];
      v[j] /= d+1;
      }
    }

  hi = new simplexHint[csiAll];
  for (auto i=0; i<csiAll; ++i)
    if (!precomputeBary(si[i], hi[i], viCentroid[i], qi, &qC, i>=csi))
      return false;
  pC = eval(qC);
  return true;
}

void terminate()
{
  delete [] pi;
  delete [] qi;
  delete [] si;
  delete [] hi;
}

// Return which simplex in si[] contains q.
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
  // Compute q's bary-coords w.r.t. each simplex in si[].
  // If one has coordinates all nonnegative, return that one.
  for (auto i=0; i<csi; ++i)
    if (computeBary(hi[i], q, w))
      return i;
  // q wasn't in any simplex, so try the ray-simplices.
  return searchRaySimplices(q, w);
}

int searchEdahiro(const vertex& q, double* w)
{
  if (d != 2) {
    printf("internal error: edahiro needs d==2.\n");
    exit(1);
  }
  const auto i = Edahiro_RegionFromPoint(q[0], q[1]);
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
  // q wasn't in any simplex, so try the ray-simplices.
  return searchRaySimplices(q, w);
}

int findSimplex(const vertex& q, double* w) {
#ifdef TESTING
  const int iS = searchBruteForce(q, w);
  if (d == 2 && searchEdahiro(q, w) != iS)
    printf("searchEdahiro and searchBruteForce disagree\n");
  return iS;
#else
  return d==2 ? searchEdahiro(q, w) : searchBruteForce(q, w);
#endif
}

#define DATAVIZ
#ifdef DATAVIZ
int iFound = 0;
double wFound[d+1] = {0};
vertex rFound{0};
#endif

// Map a d-vertex to an e-vertex.
e_vertex eval(const vertex& q)
{
  // Find which simplex s contains q.

  double w[d+1]; // q's coordinates w_j with respect to s.
  const int iS = findSimplex(q, w);
  const simplex& s = si[iS];

#ifdef TESTING
  // Verify that q == the point whose barycoords are w[] wrt s.
  vertex r{0};
  for (auto j=0; j<d; ++j) {
    for (auto i=0; i<d+1; ++i)
      r[j] += w[i] * (s[i] < 0 ? qC : qi[s[i]])[j];
  }
  // The reconstructed point is r.  How far is it from q?
  // (How accurate were the barycoords w[]?)
  auto dist = 0.0;
  for (auto j=0; j<d; ++j)
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
  for (auto j=0; j<d+1; ++j)
    wFound[j] = w[j];
#endif

  // Sum with weights w[] and vertices pi[s[]].
  e_vertex p{0};
  for (auto j=0; j<e; ++j)
    for (auto i=0; i<d+1; ++i)
      p[j] += w[i] * (s[i] < 0 ? pC : pi[s[i]])[j];
  return p;
}

#ifdef DATAVIZ

constexpr auto NaN = std::numeric_limits<double>::signaling_NaN();
vertex vQ{NaN, NaN};
e_vertex vP{0};
constexpr auto margin = 0.15 * scale;

void drawChar(const vertex& v, char c) {
  glRasterPos2f(v[0], v[1]);
  glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c);
}

// Convert std::array to GLdouble*.
void glVert2(const vertex& v) {
  glVertex2d(v[0], v[1]);
}

void drawSimplices(bool fInside) {
  if (fInside)
    glColor3f(1,0,.4);
  else
    glColor3f(.3,0,.12);
  for (auto i=0; i<csi; ++i) {
    const auto& s = si[i];
    glBegin(GL_LINE_LOOP);
    for (auto j=0; j<=d; ++j)
      glVert2(qi[s[j]]);
    glEnd();
  }
}

void drawRaysimplices(bool fInside) {
  if (fInside)
    glColor3f(0,.35,0);
  else
    glColor3f(0,0.7,0);
  for (auto i=csi; i<csiAll; ++i) {
    const auto& s = si[i];
    glBegin(GL_LINE_LOOP);
      glVert2(qC);
      glVert2(qi[s[0]]);
      glVert2(qi[s[1]]);
    glEnd();
  }
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT);

  // Mouse moved, thus assigning to vQ.
  // GLUT can't report mouse position until then.
  const bool fInited = !std::isnan(vQ[0]);

  const bool fInside = iFound < csi;
  const simplex& s = si[iFound];
  // Draw layers, most to least hidden.

  if (fInited) {
    // Bounding simplex (i.e., triangle).
    glBegin(GL_TRIANGLES);
    if (fInside) {
      glColor3f(.25,0,.08);
      for (auto j=0; j<=d; ++j)
	glVert2(qi[s[j]]);
    } else {
      glColor3f(0,.15,0);
      glVert2(qC);
      auto v = qi[s[0]];
      v[0] += 100.0 * (v[0] - qC[0]);
      v[1] += 100.0 * (v[1] - qC[1]);
      glVert2(v);
      v = qi[s[1]];
      v[0] += 100.0 * (v[0] - qC[0]);
      v[1] += 100.0 * (v[1] - qC[1]);
      glVert2(v);
    }
    glEnd();
  }

  // Bar graph of output values.
  glBegin(GL_QUADS);
  for (auto i=0; i<e; ++i) {
    constexpr auto y0 = -0.5 * margin;
    const auto y1 = vP[i];
    const auto x0 = scale / e * i;
    const auto x1 = scale / e * (i+0.6);
    glColor3f(0.0, 0.0, 0.1); glVertex2d(x0, y0);
    glColor3f(0.3, 0.3, 0.7); glVertex2d(x0, y1);
    glColor3f(0.3, 0.3, 0.7); glVertex2d(x1, y1);
    glColor3f(0.0, 0.0, 0.1); glVertex2d(x1, y0);
  }
  glEnd();

  if (fInited) {
    // Disc around each vertex of containing simplex.
    // Disc's radius indicates weight (barycentric coordinate) of corresponding vertex.
    for (auto i=0; i<=d; ++i) {
      const auto iVertex = s[i];
      const auto& v = iVertex < 0 ? qC : qi[iVertex];
      const auto r0 = 0.07 * scale;
      const auto r  = 0.07 * scale * sqrt(fabs(wFound[i]));
      glPushMatrix();
      glTranslatef(v[0], v[1], 0);
      glColor3f(0.5,0.5,0.5);
      glScalef(r, r, 0);
      // This isn't worth precomputing or moving to a display list.
      glBegin(GL_POLYGON);
	for (auto a = 0.0; a <= 2.0*M_PI; a += 0.05)
	  glVertex2f(cos(a), sin(a));
      glEnd();
      glColor3f(0.2,0.2,0.2);
      glScalef(r0/r, r0/r, 0);
      glBegin(GL_LINE_LOOP);
	for (auto a = 0.0; a <= 2.0*M_PI; a += 0.05)
	  glVertex2f(cos(a), sin(a));
      glEnd();
      glPopMatrix();
    }
  }

  // Color the hull correctly.
  if (fInside) {
    drawRaysimplices(fInside);
    drawSimplices(fInside);
  } else {
    drawSimplices(fInside);
    drawRaysimplices(fInside);
  }

  // Vertices.
  glColor3f(1,1,0);
  for (auto i=0; i<cPoint; ++i)
    drawChar(qi[i], '0' + i);

  glColor3f(1,1,1);
  drawChar(qC, 'C');

  if (fInited) {
#ifdef TESTING
    // Reconstructed point.
    glColor3f(1,.2,.2);
    drawChar(rFound, 'R');
#endif
    // Query point.
    glColor3f(1,1,1);
    drawChar(vQ, 'q');
  }

  glutSwapBuffers();
}

int xSize = 700;
int ySize = 700;

inline void XYFromMouse(double& x, double& y, int xM, int yM)
{
  x =       double(xM) / xSize  * (scale + 2*margin) - margin;
  y = (1. - double(yM) / ySize) * (scale + 2*margin) - margin;
}

void mouse_hover(int x, int y)
{
  XYFromMouse(vQ[0], vQ[1], x, y);
  vP = eval(vQ);
}

#ifdef MUST_HOLD_DOWN_MOUSE_BUTTON
void mouse(int button, int state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    mouse_hover(x, y);
}
#endif

void keyboard(unsigned char key, int /*x*/, int /*y*/)
{
  switch (key)
    {
  case 'q':
  case 27: /* escape */
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
  gluOrtho2D(-margin, margin+scale, -margin, margin+scale);
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
#ifdef MUST_HOLD_DOWN_MOUSE_BUTTON
  glutMouseFunc(mouse);
#else
  glutPassiveMotionFunc(mouse_hover);
#endif
  glutMotionFunc(mouse_hover);
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutIdleFunc(display);
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();
  glutMainLoop(); // never returns
}

#else

void evalAutomatic()
{
  // Exercise the evaluation function eval().
  // Compiled -O2 or -O3 for a 1GHz Pentium III,
  // this does one test in 310 usec or 3200 tests per second for e=42, d=7.
  // For d=3, 100 usec/test or 10000 tests/sec.
  // For d=2, 66 usec/test or 15000 tests/sec.

  constexpr auto ctest = 20; // 100000 for timing tests
  vertex qtest[ctest];
  randomSites_d(qtest, ctest);
  for (int i=0; i<ctest; ++i)
    {
    printf("eval %d/%d\n", i, ctest);
    const auto q = qtest[i];
    dump_d("query: ", q);
    const auto p = eval(q);
    dump_e("result: ", p);
    }
}
#endif

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
  if (!init())
    return -1;
#ifdef DATAVIZ
  evalInteractive(argc, argv);
#else
  evalAutomatic();
#endif
  terminate();
  return 0;
}
