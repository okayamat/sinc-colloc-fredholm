#include "matrixvector.h"

double* AllocVec(int dim)
{
  double* vec;
  vec = (double*)calloc(dim, sizeof(double));

  /* Check memory */
  if (vec == NULL) {
    free(vec);
    fprintf(stderr, "Vector was not allocated!\n");
    exit(EXIT_FAILURE);
  }

  return vec; /* Return the first address of the allocated memory */
}

void FreeVec(double* v)
{
  free(v);
}

int lapack_linsolve(double* A, double* b, int n) {
  int nrhs, lda, ldb, info;
  int* ipiv;

  nrhs = 1;
  lda = n;
  ldb = n;

  ipiv = (int*)malloc(n*sizeof(int));

  /* Check memory */
  if (ipiv == NULL) {
    free(ipiv);
    fprintf(stderr, "ipiv was not allocated!\n");
    exit(EXIT_FAILURE);
  }

  dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);

  free(ipiv);

  return info;
}
