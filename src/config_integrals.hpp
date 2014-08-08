#ifndef __CINTEGRAL_HPP_
#define __CINTEGRAL_HPP_

#include <iostream>
#include "utils.hpp"
#include "integrator.hpp"
#include "potentials.hpp"

#define INT_ORDER3 20
#define INT_ORDER6 10
//#define MAX_DIST(rho) (LJ_CUTOFF*LJ_SIGMA+0.5*ASD(rho))

using namespace std;
using namespace std::placeholders;

/*
 * function to return the heuristically chosen bound
 * for the integration of the potential
 * actual integration domain is wigner-seitz cell, but domain
 * is taken up to range of relevant contributions
 * args:
 *   po -> pointer to an instance of some derivative of the Potential class
 *   contributors -> vector of 3-length int arrays representing particles
 *                   that contribute to the potential
 *   sim -> struct of relevant simulation information
 * returns:
 *   integration bound for cubic integration domain
 */
double get_bounds(Potential *po, vector<int*> contributors, sim_i *sim);

/*
 * function to return a modified integration bound for the case of two
 * moving particles who are nearest neighbors by moving one away from the
 * other by the amount of the original bound, and using the same heuristic
 * as in the one particle or non-nn case
 * args:
 *   po -> pointer to an instance of some derivative of the Potential class
 *   cs1 -> vector of 3-length int arrays representing particles
 *          that contribute to the potential for the 1st particle
 *   cs2 -> same thing for second particle
 *   from2 -> 3-length int array containing index of second particle
 *            NOTE: if 1st particle is (l,m,n), second should be (l+1,m,n)
 *   sim -> struct of relevant simulation information
 * returns:
 *   extended integration bound for the +x direction for 1st particle
 *   and -x direction for the 2nd particle
 */
double get_nn_bound(Potential *po, vector<int*> cs1, vector<int*> cs2,
                    int *from2, sim_i *sim);

/*
 * computes the bounds as in get_bounds, but does it in every direction
 * individually (best when vacancies are around)
 * args:
 *   po -> pointer to an instance of some derivative of the Potential class
 *   cs -> vector of 3-length int array representing particles that contribute
 *         to the potential
 *   sim -> struct of relevant simulation information
 * returns:
 *   array of all integration bounds for particle. ordered by
 *   [+x,+y,+z,-x,-y,-z]
 */
double *get_all_bounds(Potential *po, vector<int*> cs, sim_i *sim);

/*
 * computes the 1st order (i.e. one moving particle in an otherwise frozen
 * matrix) configuration integral. essentially the 3d integral of the Boltzmann
 * factor over an approximation of the ws-cell
 * args:
 *   po -> pointer to ani nstance of some derivative of the Potential class
 *   sim -> struct of relevant simulation information
 *   vac1 -> optional. if given, expect 3-length array of indices representing
 *           a vacancy in the lattice
 *   vac2 -> optional. if given, expect 3-length array of indices representing
 *           a vacancy in the lattice
 * returns:
 *   value of 1st order configuration integral
 */
double config_integral_1p(Potential *po, sim_i *sim, 
                          int *vac1 = NULL, int *vac2 = NULL);

/*
 * computes the 2nd order (i.e. two moving particles in an otherwise frozen
 * matrix) configuration integral.  eseentially the 6d integral of the product
 * of the Boltzmann factors of the two particles of both ws cells
 * args:
 *   po -> pointer to an instance of some derivative of the Potential class
 *   sim -> struct of relevant simulation information
 *   from2 -> indices of second moving particle
 *   shell -> shell of the second moving particle relative to the first
 *            really only matters if nn (i.e. shell = 1).  all other
 *            values have same behavior.  
 *            NOTE: if nn, expect (l+1,m,n) if 1st particle is (l,m,n)
 * returns:
 *   value of 2nd order configuration integral
 */
double config_integral_2p(Potential *po, sim_i *sim, int *from2, int shell,
                          int *vac1 = NULL);

#endif
