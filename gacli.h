// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

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
