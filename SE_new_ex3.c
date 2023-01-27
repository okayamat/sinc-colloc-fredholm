#include "matrixvector.h"
#include "SE_basis_func.h"
#include <time.h>

double beta(double p, double q)
{
  return tgamma(p)*tgamma(q)/tgamma(p+q);
}

double a_l(int l)
{
  return pow(3.0*M_1_PI, l);
}

double b_l(int l)
{
  return pow(2*M_SQRT2/3.0, l);
}

/* kSE(x, a, b, tau) = k(x, SE_trans(a, b, tau)) */
double kSE(double x, double a, double b, double tau)
{
  double val = 0;
  double tma = (b - a) * SE_trans(0, 1, tau); /* t - a */
  double bmt = (b - a) * SE_trans(1, 0, tau); /* b - t */
  for (int l = 1; l <= 100; l++) {
    val += pow(tma, a_l(l)) * pow(bmt, 1-b_l(l));
  }

  return pow(x, sqrt(3)-1)*val;
}

double g(double x)
{
  double val = 0;
  for (int l = 1; l <= 100; l++) {
    val += beta(a_l(l)+1.5 , 2.0 - b_l(l));
  }
  val *= pow(x, sqrt(3)-1);

  return sqrt(x) - val;
}

double u(double x)
{
  return sqrt(x);
}

double uSEn(double a, double b, double x, double h, int N, double* f_N, int n)
{
  int j;
  double t = SE_trans_inv(a, b, x);
  double ans = 0;

  for (j = N; j > 0; j--) {
    ans += (f_N[ j+N] - f_N[0]*waSE( j*h) - f_N[n-1]*wbSE( j*h)) * S( j, h, t);
    ans += (f_N[-j+N] - f_N[0]*waSE(-j*h) - f_N[n-1]*wbSE(-j*h)) * S(-j, h, t);
  }
    ans += (f_N[ 0+N] - f_N[0]*waSE( 0  ) - f_N[n-1]*wbSE( 0  )) * S( 0, h, t);

    ans += f_N[0]*wa(a, b, x) + f_N[n-1]*wb(a, b, x);

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
  double* f_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    f_N[i+N] = g(x_N[i+N]);
  }

  return f_N;
}

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = 3.14;
  double alpha = 0.5;
  double *A_N, *f_N, *x_N;
  int i, n, N, info;
  double h, err, maxerr, x;
  int SAMPLE = 1000;
  double hh = (b - a)/SAMPLE;
  clock_t start, end;
  double time;

  for (N = 5; N <= 150; N += 5) {
    start = clock();

    n = 2*N+1;
    h = sqrt(M_PI*d / (alpha * N));

    x_N = substitute_xN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, x_N);
    f_N = substitute_fN(N, n, x_N);

    info = lapack_linsolve(A_N, f_N, n);

    if ( info == 0 ) {

      maxerr = 0;

      for (i = 1; i < SAMPLE; i++) {
        x = a + i*hh;

        err = fabs(u(x) - uSEn(a, b, x, h, N, f_N, n));

        maxerr = fmax(err, maxerr);
      }

      printf("%d\t%e\t", N, maxerr);
    } else fprintf(stderr, "error in lapack_linsolve!\n");

    FreeVec(f_N);
    FreeVec(A_N);
    FreeVec(x_N);

    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%e\n", time);
  }

  return EXIT_SUCCESS;
}
