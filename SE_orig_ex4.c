#include "matrixvector.h"
#include "SE_basis_func.h"

double beta(double p, double q)
{
  return tgamma(p)*tgamma(q)/tgamma(p+q);
}

/* kSE(x, a, b, tau) = k(x, SE_trans(a, b, tau)) */
double kSE(double x, double a, double b, double tau)
{
  double  t  = SE_trans(a, b, tau);
  double t18 = t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t*t; /* t^(18) */
  double tma = (b - a) * SE_trans(0, 1, tau); /* t - a */
  double bmt = (b - a) * SE_trans(1, 0, tau); /* b - t */
  double tmp = bmt * tma; /* (b - t)(t - a) */
  double val = (1+t18*t*t)*(2 + t + t18*t*t*t);
  val += (5 * (2 + x*x) * t18 * (tmp));
  val *= 2 * pow(tmp, (2 - x*x)/(2 + x*x)) / ((2+x*x)*(1+t18*t*t));

  return val;
}

double g(double x)
{
  return 2*x/(1+ pow(x, 20)) - 4*beta(1.5, 4.0/(2+x*x))/(2+x*x);
}

double u(double x)
{
  return 2*x/(1+ pow(x, 20));
}

double del0(int i, int j)
{
  if (i == j)
    return 1.0;
  else
    return 0.0;
}

double FredSEn_wa(double a, double b, double x, double h, int N)
{
  double sum=0;

  for (int j = N; j > 0; j--) {
    sum += kSE(x, a, b, j*h) * waSE( j*h) * SE_trans_div(a, b, j*h);
    sum += kSE(x, a, b,-j*h) * waSE(-j*h) * SE_trans_div(a, b,-j*h);
  }
    sum += kSE(x, a, b, 0  ) * waSE( 0  ) * SE_trans_div(a, b, 0  );

  return sum*h;
}

double FredSEn_wb(double a, double b, double x, double h, int N)
{
  double sum=0;

  for (int j = N; j > 0; j--) {
    sum += kSE(x, a, b, j*h) * wbSE( j*h) * SE_trans_div(a, b, j*h);
    sum += kSE(x, a, b,-j*h) * wbSE(-j*h) * SE_trans_div(a, b,-j*h);
  }
    sum += kSE(x, a, b, 0  ) * wbSE( 0  ) * SE_trans_div(a, b, 0  );

  return sum*h;
}

double uSEn(double a, double b, double x, double h, int N, double* f_N, int n)
{
  int j;
  double t = SE_trans_inv(a, b, x);
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
    x_N[i+N+1] = SE_trans(a, b, i*h);
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
      - FredSEn_wa(a, b, x_N[i+N+1], h, N);
  }

  for (j = -N; j<= N; j++) {
    for (i = -N-1; i <= N+1; i++) {
      A_N[i+N+1 + (j+N+1)*n] = del0(i, j)
        - h*kSE(x_N[i+N+1], a, b, j*h)*SE_trans_div(a, b, j*h);
    }
  }

  j = N+1;
  for (i = -N-1; i <= N+1; i++) {
    A_N[i+N+1 + (j+N+1)*n] = wb(a, b, x_N[i+N+1])
      - FredSEn_wb(a, b, x_N[i+N+1], h, N);
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
  double a =-1.0;
  double b = 1.0;
  double d = 1.57;
  double alpha = 1.0;
  double *A_N, *f_N, *x_N;
  int i, n, N, info;
  double h, err, maxerr, x;
  int SAMPLE = 1000;
  double hh = (b - a)/SAMPLE;

  for (N = 5; N <= 150; N += 5) {

    n = 2*N+3;
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

      printf("%d\t%e\n", N, maxerr);
    } else fprintf(stderr, "error in lapack_linsolve!\n");

    FreeVec(f_N);
    FreeVec(A_N);
    FreeVec(x_N);
  }

  return EXIT_SUCCESS;
}
