#ifndef __UTIL_HPP_
#define __UTIL_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//macro to compute the distance between nn spheres in an fcc lattice
//as a function of density.  derivation of form in my research notebook
#define ASD(rho) (cbrt(sqrt(2.)/(rho)))
//macro to compute 1/kT with Boltzmann constant as 1 for now
#define BETA(T) (1./(T))

using namespace std;

/*
 * data structure that contains info used by many functions to make
 * calling signatures easier
 */
struct sim_i{
  int *from;
  double rho;
  double asd;
  double T;
};
//typedef struct sim_info sim_i;

/*
 * function to test for equivalence between two sets of indices
 * the compiler will inline this, but I broke it out cause its used
 * in several places
 * args:
 *   p1 -> first particle index. must be length 3
 *   p2 -> second particle to compare to the first. also length 3
 * returns:
 *   true iff all entries of p1 equivalent to the corresponding ones in p2
 *   false otherwise
 */
bool index_eq(int *p1, int *p2);

/*
 * simple utility function to print out the indices of a particle
 * again the compiler will probably inline this
 * args:
 *   f -> 3 length int [] to print out
 * returns:
 *   none (yay side effects)
 */
void print_indices(int *f);

double fcc_raw(double *c, int *from, int *to, double asd, double *c2 = NULL);

/*
 * function to compute distance between two points in an fcc lattice
 * args:
 *   c -> coordinates relative to fixed lattice site for particle to compute
 *        distances from
 *   from -> 3 length int array containing particle indices to compute
 *           distances from
 *   to -> 3 length int array containing particle indices to compute distances
 *         to
 *   rho -> particle density in lattice
 *   c2 -> optional. coordinates relative to fixed lattice site for 'to' ptcl
 *         if given, distance is computed to these coordinates
 * returns:
 *   distance from c relative to 'from' to c2 (if given) relative to 'to', and
 *   to the fixed lattice site of 'to' otherwise
 */
double fcc_distance(double *c, int *from, int *to, double asd, 
                    double *c2 = NULL);

/*
 * finds all particles within a given maximum distance of the ptcl at f
 * args:
 *   f -> 3 length int[]. particle to compute distances from
 *   rho -> density of particles in lattice
 *   max_distance -> cutoff distance to find particles within
 * returns:
 *   vector of 3-length integer arrays containing the indices of all particles
 *   within max_dist of f
 */
vector<int*> find_contributors(int *f, double rho, double max_dist);

/*
 * organizes all particles within max_dist of f into groups based on distance
 * args:
 *   f -> 3 length int[]. particle indices to compute distances from
 *   rho -> density of particles in lattice
 *   max_distance -> cutoff distance to find particles within
 * returns:
 *   vector containing vectors of 3-length integer arrays.  all indices
 *   within each inner vector are at the same distance from f, sorted by
 *   distance in order of increasing index
 */
vector<vector<int*>> get_shells(int *f, double rho, double max_dist);

/*
 * overloaded version of above that allows caller to get a maximum number
 * of unique shells instead of all shells within a given distance
 * args:
 *   same as above
 * returns:
 *   same as above, except contains nm unique shells
 */
vector<vector<int*>> get_shells(int *f, double rho, int nm);

/*
 * function to find which shell a given index is in
 * args:
 *   shells -> set of indices organized into shells. c.f. return value for
 *             get_shells
 *   target -> particle indices to search for
 * returns:
 *   index into shells corresponding to vector in which target is a member
 *   -1 if there is no such index
 */
int search_shells(vector<vector<int*>> shells, int *target);

vector<int*> shell_intersect(vector<vector<int*>> shells1, 
                             vector<vector<int*>> shells2);

/*
 * function to determine whether a given point is inside the wigner-seitz
 * cell of particle with indices f
 * args:
 *   f -> 3 length int[]. particle indices to compute distances from
 *   rho -> density of particles in lattice
 *   pt -> coordinates relative to fixed lattice point of f
 * returns:
 *   true iff pt is inside the wigner-seitz cell of f
 *   false otherwise
 */
bool query_ws(int *f, double rho, double *pt);

/*
 * function to output a table containing information about particles
 * within a certain number of shells of a given particle
 * args:
 *   f -> particle indices to compute distances from
 *   rho -> density of particles in lattice
 *   nm -> maximum number of shells to consider
 * returns:
 *   none
 */
void print_ptcl_table(int *f, double rho, vector<vector<int*>> shells);

#endif
