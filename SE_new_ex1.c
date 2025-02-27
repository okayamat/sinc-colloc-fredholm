#include "matrixvector.h"
#include "SE_basis_func.h"

/* kSE(x, a, b, tau) = k(x, SE_trans(a, b, tau)) */
double kSE(double x, double a, double b, double tau)
{
  double t = SE_trans(a, b, tau);
  return x*t;
}

double g(double x)
{
  double r = 0.5;
  return r/((x-0.5)*(x-0.5) + r*r) - x*atan(1/(2*r));
}

double u(double x)
{
  double r = 0.5;
  return r/((x-0.5)*(x-0.5) + r*r);
}

double uSEn(double a, double b, double x, double h, int N, double* c_N, int n)
{
  int j;
  double t = SE_trans_inv(a, b, x);
  double ans = 0;

  for (j = N; j > 0; j--) {
    ans += c_N[ j+N] * S( j, h, t);
    ans += c_N[-j+N] * S(-j, h, t);
  }
    ans += c_N[ 0+N] * S( 0, h, t);

    ans += c_N[n]*wa(a, b, x) + c_N[n+1]*wb(a, b, x);

  return ans;
}

double* substitute_xN(double a, double b, double h, int N, int n)
{
  int i;
  double* x_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    x_N[i+N] = SE_trans(a, b, i*h);
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
        -= h*kSE(x_N[i+N], a, b, j*h)*SE_trans_div(a, b, j*h);
    }
  }

  return A_N;
}

double* substitute_fN(int N, int n, double* x_N)
{
  int i;
  double* f_N = AllocVec(n+2);

  for (i = -N; i <= N; i++) {
    f_N[i+N] = g(x_N[i+N]);
  }

  return f_N;
}

void translate_fN(double h, int N, int n, double* f_N)
{
  int j;
  f_N[n+1] = f_N[n-1];
  f_N[n]   = f_N[0];

  for (j = N; j >= -N; j--) {
    f_N[j+N] = f_N[j+N] - f_N[n]*waSE(j*h) - f_N[n+1]*wbSE(j*h);
  }
}

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = 1.57;
  double alpha = 1.0;
  double *A_N, *f_N, *x_N;
  int i, n, N, info;
  double h, err, maxerr, x;
  int SAMPLE = 1000;
  double hh = (b - a)/SAMPLE;

  for (N = 5; N <= 150; N += 5) {

    n = 2*N+1;
    h = sqrt(M_PI*d / (alpha * N));

    x_N = substitute_xN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, x_N);
    f_N = substitute_fN(N, n, x_N);

    info = lapack_linsolve(A_N, f_N, n);

    if ( info == 0 ) {

      maxerr = 0;
      translate_fN(h, N, n, f_N);

      for (i = 1; i < SAMPLE; i++) {
        x = a + i*hh;

        err = fabs(u(x) - uSEn(a, b, x, h, N, f_N, n));

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
