#include "matrixvector.h"
#include "DE_basis_func.h"

/* kDE(x, a, b, tau) = k(x, DE_trans(a, b, tau)) */
double kDE(double x, double a, double b, double tau)
{
  double t = DE_trans(a, b, tau);
  return pow(x*t, 0.75);
}

double g(double x)
{
  return sqrt(x)*(1 - M_PI*M_PI*pow(M_PI_2*x, 0.25)/9.0);
}

double u(double x)
{
  return sqrt(x);
}

double uDEn(double a, double b, double x, double h, int N, double* f_N, int n)
{
  int j;
  double t = DE_trans_inv(a, b, x);
  double ans = 0;

  for (j = N; j > 0; j--) {
    ans += (f_N[ j+N] - f_N[0]*waDE( j*h) - f_N[n-1]*wbDE( j*h)) * S( j, h, t);
    ans += (f_N[-j+N] - f_N[0]*waDE(-j*h) - f_N[n-1]*wbDE(-j*h)) * S(-j, h, t);
  }
    ans += (f_N[ 0+N] - f_N[0]*waDE( 0  ) - f_N[n-1]*wbDE( 0  )) * S( 0, h, t);

    ans += f_N[0]*wa(a, b, x) + f_N[n-1]*wb(a, b, x);

  return ans;
}

double* substitute_xN(double a, double b, double h, int N, int n)
{
  int i;
  double* x_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    x_N[i+N] = DE_trans(a, b, i*h);
  }

  return x_N;
}

double* substitute_AN(double a, double b, double h, int N, int n, double* x_N)
{
  int i, j;
  double* A_N = AllocVec(n*n); /* Column-major */

  for (i = -N; i <= N; i++) {
    A_N[i+N + (i+N)*n] = 1.0;
  }

  for (j = -N; j<= N; j++) {
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n]
        -= h*kDE(x_N[i+N], a, b, j*h)*DE_trans_div(a, b, j*h);
    }
  }

  return A_N;
}

double* substitute_fN(int N, int n, double* x_N)
{
  int i;
  double* f_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    f_N[i+N] = g(x_N[i+N]);
  }

  return f_N;
}

int main()
{
  double a = 0.0;
  double b = M_PI_2;
  double d = 1.57;
  double alpha = 0.5;
  double *A_N, *f_N, *x_N;
  int i, n, N, info;
  double h, err, maxerr, x;
  int SAMPLE = 1000;
  double hh = (b - a)/SAMPLE;

  for (N = 5; N <= 60; N += 5) {

    n = 2*N+1;
    h = log(2*d*N/alpha)/N;

    x_N = substitute_xN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, x_N);
    f_N = substitute_fN(N, n, x_N);

    info = lapack_linsolve(A_N, f_N, n);

    if ( info == 0 ) {

      maxerr = 0;

      for (i = 1; i < SAMPLE; i++) {
        x = a + i*hh;

        err = fabs(u(x) - uDEn(a, b, x, h, N, f_N, n));

        maxerr = fmax(err, maxerr);
      }

      printf("%d\t%e\n", N, maxerr);
    } else fprintf(stderr, "error in lapack_linsolve!\n");

    FreeVec(f_N);
    FreeVec(A_N);
    FreeVec(x_N);
  }

  return EXIT_SUCCESS;
}
