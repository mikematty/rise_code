#ifndef __INTEGRATOR_HPP_
#define __INTEGRATOR_HPP_

#include <math.h>
#include <iostream>
#include <functional>
#include <omp.h>

using namespace std;

/*
 * class for a Gauss-Legendre quadrature integrator
 * works for a fixed order of approximation, specified at instantiation
 * range is dynamic
 */
class gl_integrator{
  public:
    /* CONSTRUCTOR */

    /*
     * initializes an integrator that uses the Gauss-Legendre rule
     * of order n, given at instantiation
     */
    gl_integrator(int n_);

    /* DESTRUCTOR */
    ~gl_integrator();

    /* FUNCTIONS */
    
    /*
     * function that takes a pointer to a function of one variable
     * and computes the order n approximation of the integral of that function
     * over the given bounds 
     */
    double integrate(double xl, double xh, double (*fn)(double x));

    /*
     * callable wrapper for nd-dimensional integration with fixed
     * integration limits in each coordinate
     * basically the only method that will ever get called by the user
     * args:
     *   rl -> nd-length array of lower limits of integration in each var
     *   rh -> same thing for upper limits of integration
     *   fn -> function object that takes an nd-length array of doubles and
     *         returns a double
     *   nd -> number of integration dimesnions
     */
    double integrate_nd(double *rl, double *rh, 
                        function<double(double *)>fn, int nd);
    
    double* ws;
    double* xis;

  private:
    /* VARIABLES */
    int n;       //order of GL integration

    //double* ws;  //weights
    //double* xis; //abcissas

    /* FUNCTIONS */

    /*
     * generates the abcissas (roots of the nth degree Legendre polynomial)
     * and weights (c.f. Numerical Recipes) for integration using the
     * GL quadrature
     */
    void gen_weights_and_abcissas();

    /*
     * recursive helper function that actually computes the result of
     * the nd dimensional integration
     */
    double nd_integrator(double *rl, double *rh, double *rsav,
                         function<double(double *)>fn, int nd, int depth = 0);
    //double nd_integrator(double *rl, double *rh, double *rsav,
    //                     double (*fn)(double *), int nd, int depth = 0);
};

#endif
