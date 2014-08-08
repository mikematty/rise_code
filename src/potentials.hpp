#ifndef __POTENTIALS_HPP_
#define __POTENTIALS_HPP_

#include <math.h>
#include <omp.h>
#include "utils.hpp"

#define LJ_CUTOFF 6.
#define LJ_SIGMA 1.
#define LJ_EPSILON 1.

class Potential {
  public:
    /* CONSTRUCTOR */
    Potential(); //default constructor

    /* DESTRUCTOR */
    ~Potential();

    /* FUNCTIONS */
    //p.v.f.'s to be implemented by derived classes
    //see function specs below in LJ
    virtual double exptential_1p(vector<int*> cs, sim_i *sim, 
                                 double *coords) = 0;
    virtual double exptential_2p(vector<int*> cs1, vector<int*> cs2,
                                 sim_i *sim, int *from2, double *coords) = 0;

    /* VARIABLES */
    //cutoff range of the potential
    double max_dist;
};

class LJ: public Potential {
  public:
    /* CONSTRUCTOR */
    //intiialize cutoff and shifted lj potential with given cutoff
    LJ(double cutoff_);

    /* DESTRUCTOR */
    ~LJ();

    /* FUNCTIONS */

    /*
     * function to compute the cutoff lennard-jones potential between two
     * bodies separated by a distance r
     * args:
     *   rsq -> squared separation distance
     * returns:
     *   value of cutoff lj potential
     */
    double cutoff_lj(double rsq);

    /*
     * function to compute the boltzmann factor for 1 particle using the
     * cutoff lj potential
     * args:
     *   contributors -> vector of 3-length integer arrays containing particle
     *                   indices for particles that contribute to the potential
     *   sim -> struct containing relevant simulation info (c.f. utils.h)
     *   coords -> coords relative to fixed lattice site at which to compute
     * returns:
     *   value of Boltzmann factor
     */
    double exptential_1p(vector<int*> cs, sim_i *sim, double *coords);

    /*
     * function to compute product of boltzmann factors for 2 particles using
     * the cutoff lj potential
     * args:
     *   cs1 -> vector of 3-length integer arrays containing particle indices
     *          for particles that contribute to the potential for the 1st ptcl
     *   cs2 -> same thing for 2nd particle
     *   sim -> struct containing relevant simulation info (c.f. utils.h)
     *   from2 -> 3-length int array containing indices of 2nd moving particle
     *   coords -> 6-length double array containing coords relative to fixed
     *             lattice site of first particle in first three entries and
     *             of 2nd particle in the next three entries
     * returns:
     *   value of product of the Boltzmann factor of the two particles
     */
    double exptential_2p(vector<int*> cs1, vector<int*> cs2, sim_i *sim, 
                         int *from2, double *coords);
    
    /* VARIABLES */
    //double max_dist;

  private:
    /* VARIABLES */
    double cutoff, sigma, epsilon, shift;
};

//wrapper functions for exptential calls
//same specs but with additional potential argument
double wrapper_1p(vector<int*> cs, sim_i *sim, Potential *po, double *coords);

double wrapper_2p(vector<int*> cs1, vector<int*> cs2, sim_i *sim, int *from2,
                  Potential *po, double *coords);

#endif
