#include "potentials.hpp"

Potential::Potential(){};
Potential::~Potential(){};

LJ::LJ(double cutoff_):cutoff(cutoff_){
  sigma = 1.; //distance scale
  epsilon = 1.; //energy scale
  max_dist = cutoff*sigma; //maximum range of potential
  shift = -1*pow(1/cutoff, 6.);
  shift = shift*(shift+1);
}

LJ::~LJ(){};

double LJ::cutoff_lj(double rsq){
  //potential past the cutoff is 0
  if(rsq > max_dist*max_dist) return 0.;
  
  //shift potential for continuity reasons
  double attractive_term = -1*pow(sigma/rsq, 3.);
  double repulsive_term = attractive_term*attractive_term;
  return 4*epsilon*(repulsive_term+attractive_term-shift);
}

double LJ::exptential_1p(vector<int*> cs,sim_i *sim,double *coords){
  double res = 0.;
  //#pragma omp parallel for reduction(+:res)
  for(int i = 0; i < (int)cs.size(); i++)
    res += cutoff_lj(fcc_raw(coords,sim->from,cs[i],sim->asd));
  return exp(-BETA(sim->T)*res/2.);
}

double LJ::exptential_2p(vector<int*> cs1, vector<int*> cs2, sim_i *sim,
                         int *from2, double *coords){
  //don't need two separate res vars. this is left over from an earlier
  //testing stage.  it doesn't really matter though
  double res1 = 0., res2 = 0.;
  //split coords up into (x,y,z) for ptcl 1 and 2 separately
  double coords1[] = {coords[0],coords[1],coords[2]};
  double coords2[] = {coords[3],coords[4],coords[5]};

  //both ptcls should have same # of contributors b/c at same density/temp
  for(int i = 0; i < (int)cs1.size(); i++){
    //if you're at the 2nd particle, must compute distance to coords2
    if(index_eq(cs1[i],from2)) 
      res1 += cutoff_lj(fcc_raw(coords1,sim->from,cs1[i],sim->asd,coords2));
    else res1 += cutoff_lj(fcc_raw(coords1,sim->from,cs1[i],sim->asd));
    //same conditional, other direction
    if(index_eq(cs2[i],sim->from))
      res2 += cutoff_lj(fcc_raw(coords2,from2,cs2[i],sim->asd,coords1));
    else res2 += cutoff_lj(fcc_raw(coords2,from2,cs2[i],sim->asd));
  }
  return exp(-BETA(sim->T)*(res1+res2)/2.);
}

double wrapper_1p(vector<int*> cs, sim_i *sim, Potential *po, double *coords){
  //if(!query_ws(sim->from,sim->rho,coords)) return 0.0;
  //just call the exptential function for the corresponding potential
  return po->exptential_1p(cs,sim,coords);
}

double wrapper_2p(vector<int*> cs1, vector<int*> cs2, sim_i *sim, int *from2,
                  Potential *po, double *coords){
  if(!query_ws(sim->from,sim->rho,coords)) return 0.0;
  //just call the exptential function for the corresponding potential
  return po->exptential_2p(cs1,cs2,sim,from2,coords);
}
