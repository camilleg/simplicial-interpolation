// Interactive demo, using freeglut.
// https://www.opengl.org/resources/libraries/glut/spec3/spec3.html
// http://freeglut.sourceforge.net/

#include <cmath>
#include <GL/freeglut.h> // Instead of glut.h, to get glutLeaveMainLoop().

#include "si.h"

bool fInside = true;
d_simplex sFound;
vertex wFound; // Actually d+1 weights for a weighted sum.
vertex rFound;

constexpr auto NaN = std::numeric_limits<double>::signaling_NaN();
vertex vQ{NaN, NaN};
vertex vP;
double scale = NaN;
double margin = NaN;

void drawChar(const vertex& v, char c) {
  glRasterPos2f(v[0], v[1]);
  glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c);
}

// Convert std::vector to GLdouble*.
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

void display() {
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
    const auto e = vP.size();
    for (auto i=0u; i<e; ++i) {
      const auto y0 = -0.5 * margin;
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
    const auto d = vQ.size(); // Fussy.  It must be 2.
    for (auto i=0u; i<=d; ++i) {
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

void keyboard(unsigned char key, int /*x*/, int /*y*/) {
  switch (key) {
  case 'q':
  case 27: // Esc
    terminate();
#ifdef __APPLE__
    exit(0);
#else
    glutLeaveMainLoop();
#endif
  }
}

int xSize = 950;
int ySize = 950;

void XYFromMouse(double& x, double& y, int xM, int yM) {
  x =       double(xM) / xSize  * (scale + 2*margin) - margin;
  y = (1. - double(yM) / ySize) * (scale + 2*margin) - margin;
}

void mouse_hover(int x, int y) {
  XYFromMouse(vQ[0], vQ[1], x, y);
  vP = eval(vQ, &fInside, &sFound, &wFound, &rFound);
}

void reshape(int w, int h) {
  glViewport(0, 0, w, h);
  xSize = w;
  ySize = h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-margin, scale+margin, -margin, scale+margin);
}

int main(int argc, char** argv) {
  glutInit(&argc, argv); // Parse -display, -geometry 200x200+30+50, -gldebug etc.
  if (argc != 3) {
    printf("usage: %s highDim numPoints\n", argv[0]);
    return 1;
  }
  constexpr auto d = 2;
  const auto e = atoi(argv[1]);
  const auto cPoint = atoi(argv[2]);
  wFound.resize(d+1);
  if (!init(d, e, cPoint, qi_kind::spaced))
    return 1;

  // Assumes that the smallest x and smallest y are near zero.
  scale = 0.0;
  for (const auto& q: qi)
    scale = std::max(scale, std::max(q[0], q[1]));
  margin = 0.2 * scale;

#ifndef __APPLE__
  // Freeglut bug.  This and glutLeaveMainLoop() abort, because !fgState.Initialised,
  // because src/osx is empty and in particular lacks fg_init_osx.c,
  // whose fgPlatformInitialize() would have set fgState.Initialised.
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
#endif
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(xSize, ySize);
  glutCreateWindow("Simplicial Interpolator");
  glutKeyboardFunc(keyboard);
  glutPassiveMotionFunc(mouse_hover);
  glutMotionFunc(mouse_hover);
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutIdleFunc(display);
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();
  glutMainLoop();
  return 0;
}
