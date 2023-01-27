#include "SE_trans.h"

double SE_trans(double a, double b, double t)
{
/*return 0.5 * ( (b - a) * tanh(0.5*t) + (b + a) );*/
  return b / (1 + exp(-t)) + a / (1 + exp(t));
}

double SE_trans_inv(double a, double b, double x)
{
  return 2 * atanh( (2*x - (b + a))/(b - a) );
}

double SE_trans_div(double a, double b, double t)
{
  double numerator   = 0.5;
  double denominator = cosh( 0.5*t ) * cosh( 0.5*t );

  return 0.5*(b - a) * numerator / denominator;
}
