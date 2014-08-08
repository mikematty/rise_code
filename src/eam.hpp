#ifndef __EAM_HPP_
#define __EAM_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string.h>
#include "utils.hpp"
#include "potentials.hpp"
#include "cspline.hpp"

#define DMAX (numeric_limits<double>::max())
#define FtoDB 14.3997

using namespace std;

/*
 * abstract base class for various embedded atom method potentials
 * inherits from Potential abc
 */
class eam: public Potential {
  public:
    /* DEFAULT CONSTRUCTOR */
    eam();
    
    /* DESTRUCTOR */
    ~eam();
    
    /* FUNCTIONS */
    /*
     * function to compute the electronic density contributions from s e-
     * args:
     *   r -> distance from source
     * returns:
     *   value of e- density at distance r from s e-
     */
    double rhoa_s(double r);
    double rhoa_d(double r); //same for d e-

    //p.v.f.'s to be implemented separately for each eam potential
    //
    /*
     * compute electronic density at distance r
     * args:
     *   r -> separation distance
     * returns:
     *   value of e- denstiy at distance r
     */
    virtual double rho_a(double r) = 0;

    /*
     * compute pair potential for two ptcls separated by distance r
     * args:
     *   r -> separation distance
     * returns:
     *   value of pair potential at separation distance r
     */
    virtual double V(double r) = 0;

    /*
     * compute value of embedding function
     * args:
     *   rho -> electronic density
     * returns:
     *   value of embedding function with electronic denstiy rho
     */
    virtual double F(double rho) = 0;

    //see function specs for exptenital functions in potentials.hpp
    virtual double exptential_1p(vector<int*> cs, sim_i* sim, double *coords);
    virtual double exptential_2p(vector<int*> cs1, vector<int*> cs2, 
                                 sim_i *sim, int *from2, double *coords);
    
    /* VARIABLES */
    double rho0; // density at 0 pressure
    double k0, kp0; //bulk modulus and its pressure deriv. at 0 pressure

  protected:
    /* CONSTRUCTOR */
    eam(string atom_data);

    /* FUCTIONS */
    void read_atom_data(string in_name);
    void get_R_coeffs();

  private:
    /* VARIABLES */
    //material parameters for s-orbital from H.F. wave function
    vector<double> zeta_s, Cs, Rs;
    //material parameters for d-orbital from H.F. wave function
    vector<double> zeta_d, Cd, Rd;
};

/*
 * class to be used for calculating the potential energy via the 
 * embedded atom method potential presented in
 * Daw and Baskes, 1986
 */
class eamdb86: public eam {
  public:
    /* CONSTRUCTOR */
    eamdb86(string atom_data, string eam_data);

    /* DESTRUCTOR */
    ~eamdb86();

    /* FUNCTIONS */
    
    //see p.v.f. specs in base class
    virtual double rho_a(double r);
    virtual double F(double rho);
    virtual double V(double r);

  private:
    /* FUNCTIONS */
    
    /*
     * loads material parameters for this potential into instance variables
     * args:
     *   in_name -> path of input data file
     * returns:
     *   none 
     */
    void read_eam_data(string in_name);
    
    /*
     * calculates spline knots used for interpolation in potential
     * quantities based on input data
     * args/returns: none
     */
    void get_knots();

    /* VARIABLES */
    double ns, n, z0, alpha, beta;
    double a0, Esub0, B;
    int nu;
    vector<double> rho_a_knots;
    vector<double> F_knots;
    vector<double> Fpp_knots;
};

/*
 * class to be used for calculating the potential energy via the
 * embedded atom method potential presented in
 * Foiles, 1985
 */
class eamf85: public eam {
  public:
    /* CONSTRUCTOR */
    eamf85(string atom_data, string eam_data);

    /* DESTRUCTOR */
    ~eamf85();

    /* FUNCTIONS */
    //see p.v.f. function specs
    virtual double rho_a(double r);
    virtual double V(double r);
    virtual double F(double rho);

  private:
    /* FUNCTIONS */
    void read_eam_data(string in_name);

    /* VARIABLES */
    double ns, n;
    double Za1, Za2, ZRC;
    vector<double> rho_knots, F_knots, Fpp_knots;
};

/*
 * class to be used for calculating the potential energy via the
 * embedded atom method potential presedned in
 * Mishin, 1999
 */
class eamm99: public eam {
  public:
    /* CONSTRUCTOR */
    eamm99(string atom_data, string eam_data);

    /* DESTRUCTOR */
    ~eamm99();

    /* FUNCTIONS */
    //see p.v.f. function specs
    virtual double rho_a(double r);
    virtual double F(double rho);
    virtual double V(double r);

  private:
    /* FUNCTIONS */
    void read_eam_data(string in_name);

    /* VARIABLES */
    vector<double> V_knots, Vpp_knots;
    vector<double> F_knots, Fpp_knots;
    vector<double> rhoa_knots, rhoapp_knots;
    vector<double> rho_knots, r_knots;
    int n_empty;
    double V_cutoff, rhoa_cutoff;
};

#endif
