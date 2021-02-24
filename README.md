**Simplicial Interpolation**

![screenshot](./screenshot.png)

This software lets you specify corresponding points in two multidimensional spaces,
to then smoothly interpolate between the spaces.  For instance, this lets a gamepad's
two joysticks (four dimensions) control dozens or hundreds of separate parameters.

This software implements simplicial interpolation as described in the article
[Interpolated Mappings for Musical Instruments](http://camille-g.com/os02.pdf),
[Organised Sound 7(2):85‒96](http://doi.org/10.1017/S1355771802002029), © Cambridge University Press.

***Building***

On almost any Linux, or Mac OS X 10.3‒11.2,
run the interactive OpenGL demo by typing `make`.

***Running***

*Wave the mouse over the window.  If you like, click and drag even beyond the window.*

The mouse pointer `q` (for "query") is interpreted as a weighted
sum of nearby points (blue triangle).  The size of a point's gray disc shows its weight.  The center point `C` is special,
used for unbounded simplices (triangles with one edge at infinity) when `q` lies outside the convex hull
of the fixed points.

*Hit q or the escape key to exit.*

***Customizing***

Start at the bottom of [si.c++](./si.c++).
For the OpenGL demo, call `evalInteractive()`;  alternatively,
to exercise the interpolator on randomly generated data, call `evalAutomatic()`.

The dimensions of the input and output spaces are set
by the constants `d` and `e` at the top of [si.h](./si.h).

***Contact***

Camille Goudeseune
    
[http://camille-g.com](http://camille-g.com)
    
camilleg@camille-g.com

***History***

This was first published in 2002 at http://zx81.isl.uiuc.edu/interpolation/ (defunct), and revised slightly in 2009.

Copyright 2021 Camille Goudeseune,
except for included code from Ken Clarkson's [hull.shar](http://netlib.sandia.gov/voronoi/), which is copyright 1995 AT&T.
