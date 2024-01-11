## Simplicial Interpolation

![screenshot](./screenshot.png)

This software lets you specify a correspondence between points
in two Euclidean spaces ℝ<sup>*d*</sup> and ℝ<sup>*e*</sup>,
to then smoothly interpolate between them.  For instance,
this lets a mouse (*d* = 2) or a gamepad's pair of joysticks (*d* = 4)
manipulate dozens or hundreds of parameters.

It implements simplicial interpolation, as described in
[Interpolated Mappings for Musical Instruments](http://camille-g.com/os02.pdf),  
published in [Organised Sound 7(2):85‒96](http://doi.org/10.1017/S1355771802002029), © Cambridge U. Press.

This software is licensed under the [MIT License](https://mit-license.org/), © 2024 Camille Goudeseune,  
except for code from Ken Clarkson's [hull.shar](http://www.netlib.org/voronoi/), which is © 1995 AT&T
and whose license is similar to the MIT License.

### How to build and self-test

On Windows 10 or 11, install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10), using the [Ubuntu 20](https://www.microsoft.com/store/apps/9n6svws3rx71) or [Ubuntu 22](https://apps.microsoft.com/store/detail/ubuntu-2204-lts/9PN20MSR04DW?hl=en-us&gl=US) distro.

On Linux, Windows, or macOS 10.3+ (2003+), `make test`.

### How to run the demo

On Linux or Windows, `sudo apt install freeglut3-dev g++ libgl-dev libgl1-mesa-dev make`.  
On Mac, `brew install freeglut`.

-   On Windows, install and run an X server such as [GWSL](https://opticos.github.io/gwsl/) or [VcXsrv](https://sourceforge.net/projects/vcxsrv/).  
-   `make demo`, or, to specify the high dimension *e* (say, 20) and the number of points (say, 100), `./glut 20 100`  
-   Move the mouse around over the window.  If you like, zoom with the scroll wheel.  
-   To exit, hit q or the escape key.  

The mouse pointer `q` ("query") is interpreted as a weighted sum
of the corners of its surrounding triangle.
The size of each point's gray disc shows its weight.
The special center point `C` is used for an unbounded simplex
(a triangle with one edge at infinity),
when `q` lies outside the points' convex hull.

### History

This was published in 2002 at [http://zx81.isl.uiuc.edu/interpolation/](http://web.archive.org/web/20021003120921/http://zx81.isl.uiuc.edu/interpolation/) (defunct), revised slightly in 2009, and moved to GitHub in 2018.
