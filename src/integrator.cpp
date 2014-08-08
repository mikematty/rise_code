#include "integrator.hpp"

using namespace std;

gl_integrator::gl_integrator(int n_):n(n_){
  //generation leaves 0 index empty
  ws = new double[n+1];
  xis = new double[n+1];
  gen_weights_and_abcissas();
}

gl_integrator::~gl_integrator(){
  delete [] ws;
  delete [] xis;
}

double gl_integrator::integrate(double xl, double xh, double (*fn)(double x)){
  int i;
  double xr,xm,dx;
  double res = 0.;

  xm = 0.5*(xh+xl); //midpoint of integration range
  xr = 0.5*fabs(xh-xl); //width of integration range

  //compute the result using the quadrature formula
  for(i = 1; i <= n; i++){
    dx = xr*xis[i]; //scaled abcissa
    res += ws[i]*((*fn)(xm+dx));
  }

  //multiply at end as each weight must be scaled as such, and this is fewer
  //ops to just do it once at the end
  return res*xr;
}

/*
 * functions for n dimensional integration
 */

double gl_integrator::nd_integrator(double *rl, double *rh, double *rsav,
                                    function<double(double *)>fn, int nd, 
                                    int depth){
  double res = 0.;
  //get median and range for current integration variable
  double xm = 0.5*(rh[depth]+rl[depth]);
  double xr = 0.5*fabs(rh[depth]-rl[depth]);
  
  //#pragma omp parallel for reduction(+:res)
  for(int i = 1; i <= n; i++) {
    //first update the current variable being integrated over
    //#pragma omp ordered
    rsav[depth] = xm+xr*xis[i];

    //if we are doing the innermost integral, then we just evaluate the fn
    if(depth == (nd-1)) res+= ws[i]*(fn(rsav));
    //otherwise, we must evaluate the next integral
    else res += ws[i]*nd_integrator(rl, rh, rsav, fn, nd, depth+1);
  }
  return res*xr;
}

/*
 * wrapper part of recursive helper function
 */
double gl_integrator::integrate_nd(double *rl, double *rh, 
                                   function<double(double *)>fn, int nd) {
  double *rsav = new double[nd];//initialize coordinate array
  double res = nd_integrator(rl, rh, rsav, fn, nd);//compute the result
  delete [] rsav; //clean up after yourself
  return res;
}

/* 
 * note that this function assumes that the range is [-1,1], and scaling is
 * done during integration
 */
void gl_integrator::gen_weights_and_abcissas(){
  double z1,z,pp,p3,p2,p1;

  double eps = 3.0e-11; //relative precision
  int m = (n+1)/2; //basically ceil(order/2)
  
  //since the roots and weights are symmetric we onlly need to generate half
  for(int i = 1; i <= m; i++) {
    //z is an approximation for the ith root
    z = cos(M_PI*(i-0.25)/(n+0.5));

    //now use Newton's method to refine the approximation to the root
    do { //I hate the fact that I'm using a do while loop right now
      p1 = 1.0;
      p2 = 0.0;
      for(int j = 1; j <= n; j++) {
        //use recurrence relation to get Legendre ploynomial @ z
        p3 = p2;
        p2 = p1;
        p1 = ((2.*j-1.)*z*p2-(j-1.)*p3)/j;
      }
      //p1 is now the correct Legendre polynomial
      //we now compute its derivative, pp, using a known formula involving
      //p2, which is now the polynomial of degree n-1
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp; //Newton's method step
    } while (fabs(z-z1) > eps);

    xis[i] = -z; //put in the left root
    xis[n+1-i] = z; //put in its symmetric counterpart
    ws[i] = 2./((1.-z*z)*pp*pp); //compute the weight
    ws[n+1-i] = ws[i]; //and its symmetric counterpart
  }
  
  return;
}
