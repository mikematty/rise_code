#ifndef __USERVARS_HPP_
#define __USERVARS_HPP_

/************************* Temperature ***************************************/

//Temperature. Should be a positive floating point number
#define T 1.0

//For simulations with a range of temperatures, these are
//the lower bound, upper bound, and # of points, respectively
//first two should be positive doubles, last should be positive integer
#define T_RANGE_LOW 0.1
#define T_RANGE_HI  1.0
#define T_NPTS      3

/**************************** Density ****************************************/

//Density. Should be a positive floating point number
#define RHO 1.0

//For simulations with a range of densities, these
//are the lower bound, upper bound an # of points, respectively
//allowed values are same as for trange variables
#define RHO_RANGE_LOW  0.1
#define RHO_RANGE_HI   1.0
#define RHO_NPTS       10

/****************************** 3d Integration *******************************/

//order of integration for 3d (i.e. 1 particle) integrals
//should be an integer greater than 1
#define INT_ORDER3 40

//calculates bounds in all directions individually for 1p integrals if 1,
//does nothing if set to 0
//should only be 1 or 0
#define ALL_BOUNDS_1P 0

/**************************** 6d Integration *********************************/

//order of integration for 6d (i.e. 2 particle) integrals
//should be an integer greater than 1
#define INT_ORDER6 30

//in nearest neighbor case, calculate extended bound in direction of nn
//with nn moved to original bound away if 1, does nothing if 0
//should only be 1 or 0
#define NN_BOUNDS_2P 0

//calculates bounds in all directions individually for each particle
//for 2p integrals if 1, does nothing if set to 0
//should only be 1 or 0
#define ALL_BOUNDS_2P 0

/*************************** Miscellaneous ***********************************/

//throws out points outside the ws cell during integration if 1
//does nothing if set to 0
//should only be 1 or 0
#define WS_CHECK 1

//maximum number of layers of influence
//will be used if the simulation requires an nm argument (ie most 2p things)
//should be a positive integer
#define NM 1

//maximum separation between two vacancies
//useful only for vacancy vs. distance simulation.
//should be a positive integer
#define MAX_SEP 1

//allows running of the executable as ./run T rho, i.e. with temp and density
//command line args instead of using those defined above if 1, 0 does nothing
//should only be 1 or 0
#define CMDLINE_TRHO 0

#endif
