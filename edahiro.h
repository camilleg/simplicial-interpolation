// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

#include "si.h"

extern bool Edahiro_Init(const int cpt, const vertex* qi, const int csi, const simplex* si);
extern int Edahiro_RegionFromPoint(double x, double y);
