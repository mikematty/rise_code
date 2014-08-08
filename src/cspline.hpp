#ifndef __CSPLINE__
#define __CSPLINE__

#include <limits>
#include <vector>

using namespace std;

//These functions are based on those found in numerical recipes, 
//with slight modifications

/*
 * function to compute the second derivative of a function y at each
 * of the n points (xi,yi). only needs to be called once for a given
 * set of points
 *
 * args:
 *   xis -> vector of abcissas at which y(x) has been evaluated
 *          NOTE: xis must be distince, must be length at least 2
 *   yis -> vector of values of y at each xi
 *          NOTE: must be same length as xis
 *   yp1 -> value of y'(x) at xis[0]. if = DMAX, assumes y'' 0 @ x[0]
 *   ypn -> value of y'(x) at xis[n-1]. if = DMAX, assumes y'' 0 @ x[n-1]
 *
 * returns:
 *   vector containing y''(xi) values for each xi
 */
vector<double> spline_d2(vector<double> xis, vector<double> yis, 
                         double yp1, double ypn);

/*
 * function to compute an interpolated value y(x) given n evaluated points
 * yi = y(xi) and ypp = y''(xi) using cubic spline interpolation, hence
 * SPLine INTerpolation
 *
 * args:
 *   xis -> vector of abcissas at which y(x) has been evaluated.
 *          NOTE: xis must be distinct
 *   yis -> vector of values of y at each xi
 *   ypp -> vector of values of the second derivative of y at each xi
 *   x   -> point at which to evaluate y(x)
 *
 * returns:
 *   interpolated value of y(x) at the input x
 */
double splint(vector<double> xis, vector<double> yis, vector<double> ypp,
              double x);

#endif
