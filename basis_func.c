#include <math.h>

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
