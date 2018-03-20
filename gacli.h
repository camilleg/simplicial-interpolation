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

typedef struct
{
  union
    {
    short rgl[1];	// placeholder for bigger array
    double rgz[1];
    };
} Member;
const short sHuge = 0x7fff;

extern Member* GADistanceMatrix(int cptArg, int cdimSrcArg, int cdimDstArg,
	double* rgzSrc);
