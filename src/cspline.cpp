#include "cspline.hpp"

vector<double> spline_d2(vector<double> xis, vector<double> yis, 
                         double yp1, double ypn) {
  double p, qn, sig, un, *u;
  int n = (int)xis.size();
  vector<double> ypp(n);

  u = new double[n];

  // if max, then set 'natural' (0 second deriv) b.c. for x1
  if(yp1 == numeric_limits<double>::max()) ypp[0] = u[0] = 0.0;
  // otherwise use specified value of first derivative
  else {
    ypp[0] = -0.5;
    u[0] = (3./(xis[1]-xis[0]))*((yis[1]-yis[0])/(xis[1]-xis[0])-yp1);
  }

  // decomposition loop of tridiagonal algorithm. ypp and u are used
  // for temporary storage of decomposed factors
  for(int i = 1; i < n-1; i++) {
    sig = (xis[i]-xis[i-1])/(xis[i+1]-xis[i-1]);
    p = sig*ypp[i-1]+2.0;
    ypp[i] = (sig-1.0)/p;
    u[i] = (yis[i+1]-yis[i])/(xis[i+1]-xis[i]) - 
           (yis[i]-yis[i-1])/(xis[i]-xis[i-1]);
    u[i] = (6.0*u[i]/(xis[i+1]-xis[i-1])-sig*u[i-1])/p;
  }

  // check derivative on other side of evaluated range
  if(ypn == numeric_limits<double>::max()) qn = un = 0.0;
  else{
    qn = 0.5;
    un = (3.0/(xis[n-1]-xis[n-2])) * 
         (ypn-(yis[n-1]-yis[n-2])/(xis[n-1]-xis[n-2]));
  }
  ypp[n-1] = (un-qn*u[n-2])/(qn*ypp[n-2]+1.0);

  // perform backsubstitution loop of the tridiagonal algorithm
  for(int i = n-2; i >= 0; i--) ypp[i] = ypp[i]*ypp[i+1]+u[i];

  delete[] u;
  return ypp;
}

double splint(vector<double> xis, vector<double> yis, vector<double> ypp, 
              double x){
  int klo, khi, k;
  double h, b, a;

  // find correct value using method of bisection.
  // this method works optimally when the xis are at random points
  // klo and khi are bounds on the appropriate index of x
  klo = 0;
  khi = (int)xis.size()-1;
  while(khi-klo > 1) {
    k = (khi+klo)/2;
    if(xis[k] > x) khi = k;
    else klo = k;
  }

  h = xis[khi]-xis[klo];
  a = (xis[khi]-x)/h;
  b = (x-xis[klo])/h;

  // now evaluate cubic spline polynomial
  return a*yis[klo]+b*yis[khi]+((a*a*a-a)*ypp[klo]+(b*b*b-b)*ypp[khi])*(h*h)/6;
}
