// Compact "triangular" representation of a symmetric cpt-by-cpt matrix.
// Instead of foo[i][j] in rgz[cpt][cpt],
// foo[TriIJ(i,j,cpt)] in rgz[triangularNumber(cpt-1)].

inline int triangularNumber(int _)
  { return _ * (_+1) / 2; }

inline int TriIJ(int i, int j, int cpt)
  { return triangularNumber(cpt-2-i) + cpt-1-j; }
