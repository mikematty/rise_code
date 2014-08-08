#include <iostream>
#include <math.h>
#include <functional>
#include <time.h>
#include "../../src/integrator.hpp"

using namespace std;
using namespace std::placeholders;

double parabola(double x) {
  return x*x;
}

double cosine(double x) {
  return cos(x);
}

double boltzmann(double x) {
  return exp(-x);
}

double lj(double x) {
  return exp(-(pow(x,-12.)-pow(x,-6.)));
}

double spheroid(double scale, double *r){
  double x = r[0];
  double y = r[1];
  double z = r[2];
  return scale*(x*x+y*y+z*z);
}

double multi_lj(double*r){
  double dist = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  return exp(-(pow(dist,-12.)-pow(dist,-6.)));
}

double func(double *r){
  double x = r[0];
  double y = r[1];
  double z = r[2];
  return sin(x)*y+exp(cos(z));
}

int main() {
  double res;

  clock_t t;
  t = clock();
  gl_integrator integrator(20);
  clock_t runtime = clock() - t;
  cout << "initializing integrator took " << ((double)runtime)
       << " seconds" << endl;;

  //parabola test
  //res = integrator.integrate(0,2,&parabola);
  
  //cosine test
  //res = integrator.integrate(&cosine);

  //exponentiation test
  //res = integrator.integrate(&boltzmann);
  
  //lennard-jones test
  //res = integrator.integrate(&lj);

  //multidimensiional test
  double *rl = new double[3];
  rl[0] = -1;
  rl[1] = -1;
  rl[2] = -1;

  double *rh = new double[3];
  rh[0] = 2;
  rh[1] = 2;
  rh[2] = 2;

  double (*fnn)(double, double *) = &spheroid;
  auto fn = bind(fnn,2.,_1);
  
  t = clock();
  res = integrator.integrate_nd(rl,rh, fn, 3);
  runtime = clock()-t;
  cout << "integration took " << ((double)runtime) << " seconds" << endl;
  delete [] rl;
  delete [] rh;

  cout << "result of integration is " << res << endl; 

  return 0;
}
