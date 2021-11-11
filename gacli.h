typedef struct
{
  union
    {
    short rgl[1];	// placeholder for bigger array
    double rgz[1];
    };
} Member;
const short sHuge = 0x7fff;

extern Member* GADistanceMatrix(int cptArg, int cdimSrcArg, int cdimDstArg, double* rgzSrc);
