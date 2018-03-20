**Simplicial Interpolation**

![screenshot](./screenshot.png)

This software lets you specify corresponding points in two multidimensional spaces,
to then smoothly interpolate between the spaces.  For instance, this lets a gamepad's
two joysticks (four dimensions) control dozens or hundreds of separate parameters.

This software implements simplicial interpolation as described in the article
[Interpolated Mappings for Musical Instruments](http://zx81.isl.uiuc.edu/camilleg/os02.pdf),
[Organised Sound 7(2):85‒96](http://doi.org/10.1017/S1355771802002029), © Cambridge University Press.

On Linux, run the interactive OpenGL demo by typing `make`.
This also compiles and runs on OS X 10.3.9.

To tweak, start at the bottom of [si.c++](./si.c++).
For the OpenGL demo, call `evalInteractive()`;  alternatively,
to exercise the interpolator on randomly generated data, call `evalAutomatic()`.

The dimensions of the domain (input) and range (output) spaces are set
by the constants `d` and `e` at the top of [si.h](./si.h).

Contact:
    Camille Goudeseune
    
    http://camilleg-g.com/
    
    cog@illinois.edu

This was first published in 2002 at http://zx81.isl.uiuc.edu/interpolation/, and revised slightly in 2009.

Copyright 2018 Camille Goudeseune,
except for included code from Ken Clarkson's [hull.shar](http://netlib.sandia.gov/voronoi/), which is copyright 1995 AT+T.
