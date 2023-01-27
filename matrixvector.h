#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>

#if !defined(___matrixvector___)
#define ___matrixvector___

double* AllocVec(int dim);

void FreeVec(double* v);

int lapack_linsolve(double* A, double* b, int n);

#endif // ___matrixvector___
