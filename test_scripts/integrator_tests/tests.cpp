#include <math.h>

/*
 * integral for x = -1..2, y = -1..2, z = -1..2 should be 81
 */
double f1(double x, double y, double z){
  return x*x+y*y+z*z;
}

/*
 * integral for x = -.5 .. .25, y = -.5 .. .25, z = -.5 .. .25 I got 3.14416e-07
 * don't know what answer should be
 */
double f2(double x, double y, double z){
  double r = sqrt(x*x+y*y+z*z);
  return exp(-(pow(r,-12)-pow(r,-6)));
}

/*
 * integral for x = -.5 .. .5, y = -1 .. 2, z = 0 .. 5 should be 15.5503
 */
double f3(double x, double y, double z){
  return sin(x)*y+exp(cos(z));
}
