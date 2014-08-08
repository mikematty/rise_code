#ifndef __THERM_HPP_
#define __THERM_HPP_

#include <math.h>
#include <omp.h>
#include "utils.hpp"
#include "config_integrals.hpp"

#define DRHO 1e-3

/*
 * computes the Helmholz free energy per particle to 1st order (i.e. one
 * moving particle in an otherwise frozen matrix) for no vacancies
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information (c.f. utils.hpp)
 * returns:
 *   value of free energy per particle
 */
double free_energy_1p(Potential *po, sim_i *sim);

/*
 * computes the Helmholz free energy per particle to 2nd order (i.e. two
 * moving particles in an otherwise frozen matrix) for no vacancies
 * derivation of functional form in my research notebook
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information (c.f. utils.hpp)
 *   nm -> maximum number of layers of separation between two moving particles
 * returns:
 *   value of free energy per particle
 */
double free_energy_2p(Potential *po, sim_i *sim, int nm);

/*
 * uses relation p = (F/N)' * (rho^2) to compute the pressure to 1st order
 * (i.e. one moving particle in an otherwise frozen matrix) for no vacancies.
 * Here the derivative is computed with respect to density
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation informatino (c.f. utils.hpp)
 * returns:
 *   value of pressure
 */
double pressure_1p(Potential *po, sim_i *sim);

/*
 * same as 1p, except computes to 2nd order
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information
 *   nm -> number of layers for interaction range of particle
 */
double pressure_2p(Potential *po, sim_i *sim, int nm);

/*
 * computes the density at which pressure_1p is 0 using the method
 * of bisection (c.f. numerical recipes)
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information
 * returns:
 *   density for sublimation to 1st order (i.e. rho st p = 0)
 */
double sublimation_density(Potential *po, sim_i *sim);

/*
 * same as above except zeroes pressure_2p
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information
 *   nm -> number of layers for interaction range of particle
 * returns:
 *   density s.t. p = 0
 */
double sublimation_density_2p(Potential *po, sim_i *sim, int nm);

/*
 * computes change in Gibb's free energy associated with vacancy formation to
 * 1st order.  uses relation from paper
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information (c.f. utils.hpp)
 *   nm -> maximum number of shells for range of influence of vacancy
 * returns:
 *   value of formation energy
 */
double formation_energy_1p(Potential *po, sim_i *sim, int nm);

/*
 * computes an array of the difference in helmholtz free energy between
 * two vacancies separated by the index in the array + 1 shells and
 * two vacancies separated far enough apart that they are non-interacting
 * args:
 *   po -> pointer to an instance of some derived class of Potential
 *   sim -> struct containing relevant simulation information
 *   nm -> maximum number of shells for range of influence of vacancy
 * returns:
 *   array st the entry at index i is the average differene in free energy
 *   between a configuration with vacancies separated by i+1 shells and
 *   separated by effectively infinite shells
 */ 
double *vac_nrgvd(Potential *po, sim_i *sim, int nm);
#endif
