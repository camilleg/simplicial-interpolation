// Simplicial interpolation.

#include <algorithm>
#include <cmath>
#include <cstdlib> // exit()
#include <limits>
#include <random>

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

vertex qC{0}; // Constructed common point of the ray-simplices.
std::vector<vertex> qi;
e_vertex pC{0}; // What qC maps to.
std::vector<e_vertex> pi;
std::vector<d_simplex> si, siRay;
std::vector<simplexHint> hi, hiRay;

void dump_simplex(const char* prefix, const d_simplex& s) {
  cout << prefix << "\n";
  for (auto i: s) dump_v("\t", i<0 ? qC : qi[i]);
}

// Scale the inputs to the hull algorithm, which uses exact integer arithmetic.
// Outside [1e2, 1e7], ch.c++ suffers degeneracies and overshoots.
// (Hull does this itself too, with mult_up.)
constexpr auto scale = 1e6;

std::uniform_real_distribution<double> range(0.0, scale);
std::default_random_engine rng;
template<class T> void stuff(T& v) {
  constexpr size_t tSize[std::tuple_size<T>::value]{};
  for (auto i = 0u; i < std::size(tSize); ++i)
    v[i] = range(rng);
}
template<class T> void randomSites(std::vector<T>& vec, int n) {
  vec.resize(n);
  for (auto& v: vec) stuff(v);
}

e_vertex eval(const vertex&);

// Sort each simplex's vertices, to compare multiple runs for testing,
// and to simplify other testing.
bool d_simplex_compare(const d_simplex& a, const d_simplex& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}
void sort_output(std::vector<d_simplex>& rgs, bool fRay) {
  for (auto& s: rgs)
    std::sort(s.begin(), s.end() - (fRay ? 1 : 0));
  std::sort(rgs.begin(), rgs.end(), d_simplex_compare);
}

bool init()
{
  if (e < d) {
    printf("error: e (%d) < d (%d).\n", e, d);
    return false;
  }
  rng.seed(std::random_device{}());

  // Make output sites p_i.
  const auto cPoint = d==2 ? 20 : d+10;
  randomSites(pi, cPoint);

  // Get input sites q_i.
  enum { q_i_Manual, q_i_SammonsMapping, q_i_GeneticAlgorithm };
  const int q_i_kind = q_i_GeneticAlgorithm;
  switch (q_i_kind)
    {
  case q_i_Manual:
    randomSites(qi, cPoint);
    break;
  case q_i_SammonsMapping:
    computeSammon(qi, pi, scale);
    break;
  case q_i_GeneticAlgorithm:
    double tmp[e*cPoint];
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<e; ++j)
      tmp[i*e + j] = pi[i][j];
    Member* m = GADistanceMatrix(cPoint, e, d, tmp);
    qi.resize(cPoint);
    for (int i=0; i<cPoint; ++i)
    for (int j=0; j<d; ++j)
      qi[i][j] = scale * double(m->rgl[i*d + j]) / sHuge;
    free(m);
    break;
    }
  // randomSites() and computeSammon() make no errors from bogus input,
  // but delaunay_tri() reports degeneracies and Edahiro_Init() fails gracefully.
  // GADistanceMatrix() handles errors (from bogus input) simply by crashing.

  // Store a triangulation of the qi's in si.
  // callhull.c++ avoids tightly coupling hull.h to this file,
  // coupling Ken Clarkson's hull code to Camille's simplicial interpolation code.
  extern void delaunay_tri(std::vector<d_simplex>&, std::vector<d_simplex>&, int, int);
  delaunay_tri(si, siRay, d, qi.size());
  if (si.empty())
    return false;
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

  if (d == 2 && !Edahiro_Init(qi, si))
    return false;

  // Precompute some things to speed up eval().

  // Calculate qC, the centroid of the bounding box of all simplices.
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
    vertex vC{0};
    for (auto j=0; j<d; ++j) {
      for (auto i: s) vC[j] += qi[i][j]; // k<0 is possible only for siRay, not for si.
      vC[j] /= d + 1.0;
    }
    const auto h = precomputeBary(s, vC, qi, &qC, false);
    if (!h.s)
      return false;
    hi.push_back(h);
  }

  hiRay.clear();
  for (const auto& s: siRay) {
    // Accumulate into vC the centroid of s.  Pass that to precomputeBary().
    vertex vC{0};
    for (auto j=0; j<d; ++j) {
      for (auto i: s) vC[j] += (i < 0 ? qC : qi[i])[j];
      vC[j] /= d + 1.0;
    }
    const auto h = precomputeBary(s, vC, qi, &qC, true);
    if (!h.s)
      return false;
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

#define DATAVIZ
#ifdef DATAVIZ
bool fInside = false;
d_simplex sFound;
double wFound[d+1] = {0};
vertex rFound{0};
#endif

// Return the simplex that contains q.
// On error, return an arbitrary simplex.
// Into w[0...d+1], stuff q's barycentric coordinates w.r.t. that simplex.
const d_simplex& findSimplex(const vertex& q, double* w) {
#ifdef DATAVIZ
  fInside = true;
#endif
  if (d == 2) {
    // Edahiro's algorithm.  Fast.
    const int i = Edahiro_RegionFromPoint(q[0], q[1]);
    if (i >= 0) {
      if (i >= int(si.size())) {
	printf("internal error: edahiro returned out-of-range value\n");
	return si[0];
      }
      if (!computeBary(hi[i], q, w)) {
	printf("internal error: edahiro returned wrong simplex.\n");
	return si[0];
      }
      return si[i];
    }
  } else {
    // Brute force.
    // Compute q's bary-coords w.r.t. each simplex in si[].
    // If one has coordinates all nonnegative, return that one.
    for (const auto h: hi)
      if (computeBary(h, q, w))
	return *h.s;
  }
  // q wasn't in any simplex, so look in the ray-simplices.
#ifdef DATAVIZ
  fInside = false;
#endif
  for (const auto h: hiRay)
    if (computeBary(h, q, w, true))
      return *h.s;
  // This should be impossible, because the ray-simplices cover R^d.
  // So arbitrarily choose the first simplex.
  printf("internal error in findSimplex\n");
  (void)computeBary(hi[0], q, w);
  return siRay[0];
}

// Map a d-vertex to an e-vertex.
e_vertex eval(const vertex& q)
{
  // Find which simplex s contains q.
  double w[d+1]; // q's coordinates w_j with respect to s.
  const d_simplex& s = findSimplex(q, w);

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
  sFound = s;
  std::copy(w, w + d+1, wFound); // For discs.
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
  for (const auto& s: si) {
    glBegin(GL_LINE_LOOP);
    for (auto i: s) glVert2(qi[i]);
    glEnd();
  }
}

void drawRaysimplices(bool fInside) {
  if (fInside)
    glColor3f(0,.35,0);
  else
    glColor3f(0,0.7,0);
  for (const auto& s: siRay) {
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
  // gluOrtho2D set the range for both x and y to -margin .. scale+margin.

  // Mouse moved, thus assigning to vQ.
  // GLUT can't report mouse position until then.
  const bool fInited = !std::isnan(vQ[0]);

  // Draw layers, most to least hidden.

  if (fInited) {
    // Bounding simplex (i.e., triangle).
    glBegin(GL_TRIANGLES);
    if (fInside) {
      glColor3f(.25,0,.08);
      for (auto i: sFound) glVert2(qi[i]);
    } else {
      glColor3f(0,.15,0);
      glVert2(qC);
      auto v = qi[sFound[0]];
      v[0] += 100.0 * (v[0] - qC[0]);
      v[1] += 100.0 * (v[1] - qC[1]);
      glVert2(v);
      v = qi[sFound[1]];
      v[0] += 100.0 * (v[0] - qC[0]);
      v[1] += 100.0 * (v[1] - qC[1]);
      glVert2(v);
    }
    glEnd();

    // Bar graph of output values.
    // (If e exceeds the window's width in pixels,
    // aliasing amusingly omits some bars.)
    glBegin(GL_QUADS);
    for (auto i=0; i<e; ++i) {
      constexpr auto y0 = -0.5 * margin;
      const auto y1 = vP[i];
      const auto x0 = scale * (i+0.4)/e;
      const auto x1 = scale * (i+0.6)/e;
      glColor3f(0.0, 0.0, 0.1); glVertex2d(x0, y0);
      glColor3f(0.3, 0.3, 0.7); glVertex2d(x0, y1);
      glColor3f(0.3, 0.3, 0.7); glVertex2d(x1, y1);
      glColor3f(0.0, 0.0, 0.1); glVertex2d(x1, y0);
    }
    glEnd();

    // Disc around each vertex of containing simplex.
    // Disc's radius indicates weight (barycentric coordinate) of corresponding vertex.
    for (auto i=0; i<=d; ++i) {
      const auto iVertex = sFound[i];
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
    drawRaysimplices(true);
    drawSimplices(true);
  } else {
    drawSimplices(false);
    drawRaysimplices(false);
  }

  // Vertices.
  glColor3f(1,1,0);
  auto i=0; for (const auto& q: qi) drawChar(q, '0' + i++);

  // Center vertex of ray simplices.
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
  gluOrtho2D(-margin, scale+margin, -margin, scale+margin);
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
  std::vector<vertex> qtest;
  randomSites(qtest, ctest);
  for (const auto& q: qtest) {
    dump_v("query: ", q);
    const auto p = eval(q);
    dump_v("result: ", p);
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
