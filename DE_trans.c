#include "DE_trans.h"

double DE_trans(double a, double b, double t)
{
/*return 0.5 * ( (b - a) * tanh( M_PI_2 * sinh(t) ) + (b + a) );*/
  return b / (1 + exp(-M_PI*sinh(t))) + a / (1 + exp(M_PI*sinh(t)));
}

double DE_trans_inv(double a, double b, double x)
{
  return asinh( M_2_PI * atanh( (2*x - (b + a))/(b - a) ) );
}

double DE_trans_div(double a, double b, double t)
{
  double numerator   = M_PI_2 * cosh(t);
  double denominator = cosh( M_PI_2 * sinh(t) ) * cosh( M_PI_2 * sinh(t) );

  return 0.5*(b - a) * numerator / denominator;
}
