#include "therm_quantities.hpp"

double free_energy_1p(Potential *po, sim_i *sim){
  // just use partition function = exp(-\beta F)
  return -log(config_integral_1p(po,sim))/BETA(sim->T);
}

double free_energy_2p(Potential *po, sim_i *sim, int nm){
  double res = 0.;
  double config2p;
  //a particular first nearest neighbor used for computing special bounds
  int first_nn[] = {(sim->from)[0]+1,(sim->from)[1],(sim->from)[2]};
  
  vector<vector<int*>> shells = get_shells(sim->from,sim->rho,nm);
  
  double zs = config_integral_1p(po,sim); //1st order configuration integral
  //cout << "0:"<< setprecision(10) << -log(zs)/BETA(sim->T) << ' ';
  
  //can compute each correction to the 1p free energy in parallel
  #pragma omp parallel for private(config2p) reduction(+:res)
  for(int i = 0; i < nm; i++){
    //use special first neighbor for bound computation
    if (i == 0) config2p = config_integral_2p(po,sim,first_nn,i+1);
    else config2p = config_integral_2p(po,sim,(shells[i])[0],i+1);
   
    //cout << i+1 << ":" << setprecision(10)
    //     <<-.5*(double)(shells[i].size())*log(config2p/(zs*zs))/BETA(sim->T)
    //     << ' ';
    
    res += (double)((shells[i]).size())*log(config2p/(zs*zs));
  }
  //detailed derivation of this form is in my notebook
  res = (-log(zs) - 0.5*res)/BETA(sim->T);
  //cout << nm+1 << ":" << setprecision(10) << res << endl;

  return res;
}

double pressure_1p(Potential *po, sim_i *sim){
  //use f'(x) \approx (f(x+\delta)-f(x-\delta))/(2\delta)
  double fplus=(sim->rho+=DRHO,sim->asd=ASD(sim->rho),free_energy_1p(po,sim));
  double fminus = (sim->rho -= 2*DRHO, sim->asd = ASD(sim->rho),
                   free_energy_1p(po,sim));
  sim->rho += DRHO; //set back to original
  sim->asd = ASD(sim->rho);
  //(F/N)' * rho^2 is the pressure
  return ((fplus-fminus)/(2*DRHO))*(sim->rho)*(sim->rho); 
}

double pressure_2p(Potential *po, sim_i *sim, int nm){
  //same as above, just with a free_energy_2p call
  double fplus = (sim->rho += DRHO,sim->asd = ASD(sim->rho),
                  free_energy_2p(po,sim,nm));
  double fminus = (sim->rho -= 2*DRHO,sim->asd = ASD(sim->rho),
                   free_energy_2p(po,sim,nm));
  sim->rho += DRHO;
  sim->asd = ASD(sim->rho);
  return((fplus-fminus)/(2*DRHO))*(sim->rho)*(sim->rho);
}

double sublimation_density(Potential *po, sim_i *sim){
  double drho, rho_mid, eps = 1e-10; //tolerance for change to guess
  double rho_min = 0.6; //needs to be manually entered so p < 0 at rho_min
  double rho_max = 1.2; //needs to be entered st p > 0 at rho_max
  //pressure at lower bound
  double p = pressure_1p(po,(sim->rho = rho_min, sim->asd=ASD(sim->rho),sim));
  //pressure at upper bound
  double pmid = pressure_1p(po,(sim->rho=rho_max,sim->asd=ASD(sim->rho),sim));
  double rtb = (p < 0.) ? (drho=rho_max-rho_min,rho_min):
                          (drho=rho_min-rho_max,rho_max);
  for(int i = 1; i <= 40; i++){
    sim -> rho = (rho_mid = rtb + (drho *= 0.5));
    sim -> asd = ASD(sim->rho);
    pmid = pressure_1p(po,sim);
    if(pmid <= 0.0) rtb = rho_mid; //if pressure < 0, then increase lower bnd
    //check if found root (you won't) or changes to rho have become smaller
    //than the threshold
    if(fabs(drho) < eps || pmid == 0.0) return rtb;
  }
  //should never get here unless something went terribly wrong
  cout << "well this is going downhill fast..." << endl;
  return 0.0;
}

double sublimation_density_2p(Potential *po, sim_i *sim, int nm){
  //same as above, just with calls to pressure_2p
  double drho, rho_mid, eps = 1e-10;
  double rho_min = 0.5;
  double rho_max = 1.5;
  sim -> rho = rho_min;
  sim -> asd = ASD(sim -> rho);
  double p = pressure_2p(po,sim,nm);
  double pmid=pressure_2p(po,(sim->rho=rho_max,sim->asd=ASD(rho_max),sim),nm);
  double rtb = (p < 0.) ? (drho = rho_max-rho_min,rho_min):
                          (drho = rho_min-rho_max,rho_max);
  for(int i = 1; i <= 40; i ++){
    sim -> rho = (rho_mid = rtb+(drho *= 0.5));
    sim -> asd = ASD(sim->rho);
    pmid = pressure_2p(po,sim,nm);
    cout << "guess rho = " << sim->rho << " with p = " << pmid << endl;
    if(pmid <= 0.0) rtb = rho_mid;
    if(fabs(drho) < eps || pmid == 0.0) return rtb;
  }
  return 0.0;
}

#if 0
double sublimation_density(Potential *po, sim_i *sim){
  double rho_min = 0.001;
  double rho_max = 0.1;
  sim -> rho = (rho_min+rho_max)/2.;
  double eps = 1e-5;
  double p = 2*eps;
  cout << "T is " << sim->T << endl;
  while(fabs(p) > eps){
    cout << "guess is " << setw(20) << setprecision(20) 
         << sim->rho << " with p = " << p << endl;
    p = pressure_1p(po, sim);
    if(p > 0) rho_max = sim->rho;
    else rho_min = sim->rho;
    sim->rho = (rho_max+rho_min)/2.;
  }
  cout << endl << endl;
  return sim->rho;
}
#endif

double formation_energy_1p(Potential *po, sim_i *sim, int nm){
  double zsvi, nrg = 0.0;

  //need to pick one vacancy from each shell
  vector<vector<int*>> shells = get_shells(sim->from,sim->rho,nm);

  //now need to compute configuration integrals
  double zs = config_integral_1p(po,sim); //no vacancy integral
  for(int i = 0; i < nm ; i++){
    //maybe need to think about more creative bound calculations here
    zsvi = config_integral_1p(po,sim,shells[i][0]);
    nrg += ((double)(shells[i]).size())*log(zsvi/zs);
  }
  //detailed derivation of form can be found in notebook
  return (pressure_1p(po,sim)/(sim->rho))-(sim->T)*nrg;
}

double *vac_nrgvd(Potential *po, sim_i *sim, int nm){
  int max_sep = 1; //maximum number of layers between two vacancies
  double *res = new double[max_sep];
  double *zsvis = new double[nm];
  double dF;
  int vac1[] = {0,0,0}; //put first vacancy at origin
  vector<vector<int*>> shells2;
  vector<vector<int*>> shells1 = get_shells(vac1,sim->rho,nm);
  vector<vector<int*>> net_shells = get_shells(vac1,sim->rho,max_sep);
  vector<int*> overlap;
  
  //zero vacancy configuration integral
  double zs = config_integral_1p(po,sim);
  //single vacancy configuration integrals up to the maximum separation
  for(int i = 0;i < nm;i++) zsvis[i]=config_integral_1p(po,sim,shells1[i][0]);
  
  //go layer by layer up to the maximum separation
  for(int i = 0; i < max_sep; i++){
    dF = 0.0; //reset count
    for(int j = 0; j < (int)(net_shells[i]).size(); j++){
      //net_shells[i][j] is the 2nd vacancy
      //will look at each ptcl in each relevant shell as 2nd vac
      shells2 = get_shells(net_shells[i][j],sim->rho,nm);
      overlap = shell_intersect(shells1,shells2); //ptcls in range of both
      
      for(int l = 0; l < (int)overlap.size(); l++){
        sim->from = overlap[l];//look at each ptcl in overlap as moving ptcl
        //layer of ptcl in overlap from vac 1 and 2 respectively
        int from1 = search_shells(shells1,overlap[l]);
        int from2 = search_shells(shells2,overlap[l]);
        //calculate the integral with both vacancies
        double zsvij = config_integral_1p(po,sim,vac1,net_shells[i][j]);
        //detailed derivation of form in notebook
        dF += -log((zs*zsvij)/(zsvis[from1]*zsvis[from2]));
      }
    }
    //normalize to shell size.  we should maybe rethink this
    dF /= (double)(net_shells[i]).size();
    //dF *= sim -> T; //and there's an inverse beta in there
    res[i] = dF;
  }
  delete[] zsvis;
  return res;
}
