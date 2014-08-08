#include "../../src/cspline.hpp"
#include <iostream>

using namespace std;

double y(double x){
  return x*x;
}

double yp(double x){
  return 2*x;
}

int main(){

  int n = 20;

  double *xis = new double[n];
  double *yis = new double[n];
  double *ypp = new double[n];

  for(int i = 0; i < n; i++){
    xis[i] = (double)i/2.;
    yis[i] = y(xis[i]);
  }


  spline_d2(xis,yis,n,yp(xis[0]),yp(xis[n-1]),ypp);
  for(int i = 0; i < n; i++){
    cout << "(interpolated,actual) value at " << (double)i/2+0.25 << " is ("
         << splint(xis,yis,ypp,n,(double)i/2+0.25) << ","
         << y((double)i/2+0.25) << ")" << endl;
  }

  delete[] xis;
  delete[] yis;
  delete[] ypp;
  return 0;
}
