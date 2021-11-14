typedef struct {
  union {
    short rgl[1]; // placeholder for bigger array
    double rgz[1];
  };
} Member;
Member* GADistanceMatrix(int cptArg, int cdimSrcArg, int cdimDstArg, double* rgzPt);

const short sHuge = 0x7fff;
