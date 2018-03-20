// This library implements simplicial interpolation as described in
// "Interpolated Mappings for Musical Instruments", Organised Sound 7(2),
// Cambridge University Press.  Copyright 2002 Camille Goudeseune.

// Compute the determinant of the n*n matrix m,
// with the recursive "expansion by minors" method.
 
double determinant(const double* m, const int n)
{
  if (n==2)
    return m[0] * m[3] - m[1] * m[2];

  /* if (n==3)
    return aei + cdh + bfg - ceg - bdi - afh; where abc is first row.  */
 
  const int n1 = n-1;
  double *m1 = new double[n1 * n1]; // Storage for minor matrices.
  double acc=0;
 
  for (int j1=0; j1<n; ++j1)      // For each row of m...
    {
    // Compute m1.
    int k=n, k1=0;
    for (int i=1; i<n; ++i)       // Skip first column (i==0).
    for (int j=0; j<n; ++j,++k)
      if (j != j1)                // Skip j1'th row.
        m1[k1++] = m[k];          // k == n*i + j
 
    // Recurse on m1.
    const double t = determinant(m1, n1);
 
    // Accumulate.
    acc += m[j1] * (j1%2 ? -t : t);
    }
  delete [] m1;
  return acc;
}

/*
  TODO:
  Convert this from recursion to iteration, by building a stack
  of n1 n1*n1 matrices.  Wastes a little stack space, but avoids
  function-call overhead, uses only one new/delete, and is friendly
  to optimizing compilers.

  Would diagonalizing the matrix and multiplying its eigenvalues
  be any faster?  (Convert it to upper triangular form, avoiding
  numerical instabilities via Gauss-Jordan, then multiply the diagonal entries.)
  e.g., CLHEP 1.8.0.0 Matrix/Matrix.cc HepMatrix::dfact_matrix()
  That's O(n^3), same as what we have here.
*/
