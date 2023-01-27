#include "SE_basis_func.h"

double S(int j, double h, double x)
{
  double val = M_PI*(x - j*h)/h;
  if (val == 0)
    return 1.0;
  else
    return sin(val)/val;
}

double wa(double a, double b, double x)
{
  return (b-x)/(b-a);
}

double wb(double a, double b, double x)
{
  return (x-a)/(b-a);
}

/* waSE(x) = wa(a, b, SE_trans(a, b, x)) */
double waSE(double x)
{
  return SE_trans(1, 0, x);
}

/* wbSE(x) = wb(a, b, SE_trans(a, b, x)) */
double wbSE(double x)
{
  return SE_trans(0, 1, x);
}
