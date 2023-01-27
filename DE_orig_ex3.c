#include "matrixvector.h"
#include "DE_basis_func.h"
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

/* kDE(x, a, b, tau) = k(x, DE_trans(a, b, tau)) */
double kDE(double x, double a, double b, double tau)
{
  double val = 0;
  double tma = (b - a) * DE_trans(0, 1, tau); /* t - a */
  double bmt = (b - a) * DE_trans(1, 0, tau); /* b - t */
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

double del0(int i, int j)
{
  if (i == j)
    return 1.0;
  else
    return 0.0;
}

double FredDEn_wa(double a, double b, double x, double h, int N)
{
  double sum=0;

  for (int j = N; j > 0; j--) {
    sum += kDE(x, a, b, j*h) * waDE( j*h) * DE_trans_div(a, b, j*h);
    sum += kDE(x, a, b,-j*h) * waDE(-j*h) * DE_trans_div(a, b,-j*h);
  }
    sum += kDE(x, a, b, 0  ) * waDE( 0  ) * DE_trans_div(a, b, 0  );

  return sum*h;
}

double FredDEn_wb(double a, double b, double x, double h, int N)
{
  double sum=0;

  for (int j = N; j > 0; j--) {
    sum += kDE(x, a, b, j*h) * wbDE( j*h) * DE_trans_div(a, b, j*h);
    sum += kDE(x, a, b,-j*h) * wbDE(-j*h) * DE_trans_div(a, b,-j*h);
  }
    sum += kDE(x, a, b, 0  ) * wbDE( 0  ) * DE_trans_div(a, b, 0  );

  return sum*h;
}

double uDEn(double a, double b, double x, double h, int N, double* f_N, int n)
{
  int j;
  double t = DE_trans_inv(a, b, x);
  double ans = 0;

  for (j = N; j > 0; j--) {
    ans += f_N[ j+N+1] * S( j, h, t);
    ans += f_N[-j+N+1] * S(-j, h, t);
  }
    ans += f_N[ 0+N+1] * S( 0, h, t);

    ans += f_N[0]*wa(a, b, x) + f_N[n-1]*wb(a, b, x);

  return ans;
}

double* substitute_xN(double a, double b, double h, int N, int n)
{
  int i;
  double* x_N = AllocVec(n);

    x_N[0] = a;
  for (i = -N; i <= N; i++) {
    x_N[i+N+1] = DE_trans(a, b, i*h);
  }
    x_N[n-1] = b;

  return x_N;
}

double* substitute_AN(double a, double b, double h, int N, int n, double* x_N)
{
  int i, j;
  double* A_N = AllocVec(n*n); /* Column-major */

  j = -N-1;
  for (i = -N-1; i <= N+1; i++) {
    A_N[i+N+1 + (j+N+1)*n] = wa(a, b, x_N[i+N+1])
      - FredDEn_wa(a, b, x_N[i+N+1], h, N);
  }

  for (j = -N; j<= N; j++) {
    for (i = -N-1; i <= N+1; i++) {
      A_N[i+N+1 + (j+N+1)*n] = del0(i, j)
        - h*kDE(x_N[i+N+1], a, b, j*h)*DE_trans_div(a, b, j*h);
    }
  }

  j = N+1;
  for (i = -N-1; i <= N+1; i++) {
    A_N[i+N+1 + (j+N+1)*n] = wb(a, b, x_N[i+N+1])
      - FredDEn_wb(a, b, x_N[i+N+1], h, N);
  }

  return A_N;
}

double* substitute_fN(int N, int n, double* x_N)
{
  int i;
  double* f_N = AllocVec(n);

  for (i = -N-1; i <= N+1; i++) {
    f_N[i+N+1] = g(x_N[i+N+1]);
  }

  return f_N;
}

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = 1.57;
  double alpha = 0.5;
  double *A_N, *f_N, *x_N;
  int i, n, N, info;
  double h, err, maxerr, x;
  int SAMPLE = 1000;
  double hh = (b - a)/SAMPLE;
  clock_t start, end;
  double time;

  for (N = 5; N <= 60; N += 5) {
    start = clock();

    n = 2*N+3;
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
