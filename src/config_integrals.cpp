#include "config_integrals.hpp"

double get_bounds(Potential *po, vector<int*> contributors, sim_i *sim) {
  //start measuring from the fixed lattice site
  double* coords = new double[3]();
  
  //choose step size as one thousandth of the distance to a nn
  double stepsize = (sim->asd)/1e3;
  double prev_potential, cur_potential;
  // I don't like do/while loops.  I should fix this.
  do {
    cur_potential = po->exptential_1p(contributors,sim,coords);
    if(coords[1] < 1e-4) prev_potential = cur_potential;
    coords[1] += stepsize;
  } while(((cur_potential/prev_potential)>1e-5)&&(coords[1] < ASD(sim->rho)));
  //y component should have largest bound so use that for everything
  double res = coords[1];
  delete [] coords; //cleanup
  return res;
}

double get_nn_bound(Potential *po, vector<int*> cs1, vector<int*> cs2, 
                    int *from2, sim_i *sim){
  double bound_1p = get_bounds(po,cs1,sim); //initial guess at bound
  //coords for 2nd p moved to initial bound for both particles
  double coords[] = {0.,0.,0.,bound_1p,0.,0.};
  double coords_1p[] = {bound_1p,0.,0.}; //and for just the second
  double stepsize = (sim->asd)/1e3;
  //construct a second simulation structure for the 2nd particle
  sim_i *sim2 = new sim_i;
  sim2->T = sim->T;
  sim2->rho = sim->rho;
  sim2->asd = ASD(sim->rho);
  sim2->from = from2;
  double prev_potential, cur_potential;
  do {
    //compute potential at 1st particle w/ 2nd moved over
    cur_potential = po->exptential_2p(cs1,cs2,sim,from2,coords);
    cur_potential /= po->exptential_1p(cs2,sim2,coords_1p);
    if(coords[0] < 1e-4) prev_potential = cur_potential;
    coords[0] += stepsize;
  } while(((cur_potential/prev_potential)>1e-5));//&&(coords[0]<ASD(sim->rho)));
  double res = coords[0]; //modified x-bound
  delete sim2;
  return res;
}

double *get_all_bounds(Potential *po, vector<int*> cs, sim_i *sim){
  double stepsize = (sim->asd)/1e3;
  double *res = new double[6];
  double coords[] = {0.,0.,0.};
  double prev_potential, cur_potential;
  //each positive direction
  for(int i = 0; i < 3; i++){
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
    do{
      cur_potential = po->exptential_1p(cs,sim,coords);
      if(fabs(coords[i]) < 1e-4) prev_potential = cur_potential;
      coords[i] += stepsize;
    } while(((cur_potential/prev_potential)>1e-5)&&(coords[i]< sim->asd));
    res[i] = coords[i];
  }
  //each negative direction
  for(int i = 0; i < 3; i++){
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
    do{
      cur_potential = po->exptential_1p(cs,sim,coords);
      if(fabs(coords[i]) < 1e-4) prev_potential = cur_potential;
      coords[i] -= stepsize;
    } while(((cur_potential/prev_potential)>1e-5)&&(fabs(coords[i])<sim->asd));
    res[i+3] = coords[i];
  }
  return res;
}

double config_integral_1p(Potential *po, sim_i* sim, int* vac1, int* vac2){
  //indicator variables for vacancy index positions
  int found1 = -1;
  int found2 = -1;
  gl_integrator integrator(INT_ORDER3); //initialize an integrator

  vector<int*> cs = find_contributors(sim->from,sim->rho,po->max_dist);
  double bound = get_bounds(po,cs,sim);
  //if there is a first vacancy, delete from the contributors if its there
  if(vac1 != NULL){
    for(int i = 0; i < ((int)cs.size()); i++){
      if(index_eq(cs[i],vac1)) found1 = i;
    }
    if(found1 >= 0) cs.erase(cs.begin()+found1);
  }
  //if there is a second vacancy, delete from the contributors if its there
  if(vac2 != NULL){
    for(int j = 0; j < ((int)cs.size()); j++){
      if(index_eq(cs[j],vac2)) found2 = j;
    }
    if(found2 >= 0) cs.erase(cs.begin()+found2);
  }
  //partially apply the wrapper to the exp(-\beta \phi) function
  //to known arguments to match the integrator call signature
  double(*fnn)(vector<int*>,sim_i*,Potential*,double*) = &wrapper_1p;
  auto fn = bind(fnn, cs, sim, po, _1);
  
  //set bounds and compute integral
  //bound = cbrt(1./(2.*(sim->rho)));
  double rh[3] = {bound,bound,bound}; //upper integration bounds
  double rl[3] = {-bound,-bound,-bound}; //lower integration bounds
  //compute result
  double res = integrator.integrate_nd(rl,rh,fn,3);
  return res;
}

double config_integral_2p(Potential *po, sim_i* sim, int* from2, int shell,
                          int *vac1){
  //find contributors to each moving particle
  vector<int*> cs1=find_contributors(sim->from,sim->rho,po->max_dist);
  vector<int*> cs2=find_contributors(from2,sim->rho,po->max_dist);
  
  double bound = get_bounds(po,cs1,sim);
  //bound = cbrt(1./(2.*(sim->rho)));
  double rh[6] = {bound,bound,bound,bound,bound,bound};
  double rl[6] = {-bound,-bound,-bound,-bound,-bound,-bound};
  //if there is a vacancyk, delete it if it is a contributor
  if(vac1 != NULL){
    int found1,found2;
    for(int i = 0; i < ((int)cs1.size()); i++){
      if(index_eq(cs1[i],vac1)) found1 = i;
      if(index_eq(cs2[i],vac1)) found2 = i;
    }
    if(found1 >= 0) cs1.erase(cs1.begin()+found1);
    if(found2 >= 0) cs2.erase(cs2.begin()+found2);
#if 1
    double *bounds1 = get_all_bounds(po,cs1,sim);
    int *from = sim->from;
    sim -> from = from2;
    double *bounds2 = get_all_bounds(po,cs2,sim);
    sim -> from = from;
    //WARNING: ugly hardcoding ahead
    rh[0] = bounds1[0]; //cout << "x1max = " << bounds1[0] << endl;
    rh[1] = bounds1[1]; //cout << "y1max = " << bounds1[2] << endl;
    rh[2] = bounds1[2]; //cout << "z1max = " << rh[2] << endl;
    rh[3] = bounds2[0]; //cout << "x2max = " << rh[3] << endl;
    rh[4] = bounds2[1]; //cout << "y2max = " << rh[4] << endl;
    rh[5] = bounds2[2]; //cout << "z2max = " << rh[5] << endl;
    rl[0] = bounds1[3]; //cout << "x1min = " << rl[0] << endl;
    rl[1] = bounds1[4];
    rl[2] = bounds1[5];
    rl[3] = bounds2[3];
    rl[4] = bounds2[4];
    rl[5] = bounds2[5]; 
#endif
  }

  gl_integrator integrator(INT_ORDER6);
  //partially apply wrapper to exp(-\beta(\phi_1+\phi_2)) function to known
  //arguments to match integrator calling signature
  double(*fnn)(vector<int*>,vector<int*>,
               sim_i*,int*,Potential*, double*) = &wrapper_2p;
  auto fn = bind(fnn, cs1, cs2, sim, from2, po, _1);

  //set bounds, and modify if in first shell
#if 1
  if(shell == 1){
    double bound_nn= get_nn_bound(po,cs1,cs2,from2,sim);
    rh[0] = bound_nn;
    rl[3] = -bound_nn;
  }
#endif
  //compute result
  return integrator.integrate_nd(rl,rh,fn,6);
}
