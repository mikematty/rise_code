#include "simulate.hpp"

void write_eam_info_r(eam *po, double rmin, double rmax, double nsteps){
  ofstream outf;
  outf.open("radial_scaling_db86.txt");

  outf << setw(5) << "r" << ' '
       << setw(12) << "rho_a(r)" << ' '
       << setw(12) << "V(r)" << ' '
       << setw(12) << "F(rho_a(r))" << endl;

  for(double r = rmax; r > rmin; r -= (rmax-rmin)/nsteps){
    outf << setw(5) << r << ' '
         << setw(12) << po->rho_a(r) << ' '
         << setw(12) << po->V(r) << ' '
         << setw(12) << po->F(po->rho_a(r)) << endl;
  }
  outf.close();
  return;
}

void write_eam_info_rho(eam *po, double rhomin, double rhomax, double step){
  int np;
  vector<vector<int*>> layers;
  vector<int*> cs;
  double rho_h, phi, r=0.;
  int f[] = {0,0,0};
  double origin[] = {0.,0.,0.};
  
  ofstream outf;
  outf.open("density_scaling_f85.txt");

  outf << setw(7) << "rho" << ' '
       << setw(12) << "a" << ' '
       << setw(3) << "nl" << ' '
       << setw(4) << "np" << ' '
       << setw(16) << "rho_a @ O" << ' '
       << setw(16) << "V @ O" << ' '
       << setw(16) << "F @ O" << ' '
       << setw(16) << "E_coh" << endl;

  cout.precision(10);

  for(double rho = rhomin; rho < rhomax; rho += step){
    layers = get_shells(f, rho, po->max_dist);
    np = 0;
    for(int i = 0; i < (int)layers.size(); i++) np += (int)(layers[i]).size();

    rho_h = phi = 0.;
    cs = find_contributors(f, rho, po->max_dist);
    for(int i = 0; i < (int)cs.size(); i++){
      r = fcc_distance(origin,f,cs[i],ASD(rho));
      rho_h += po->rho_a(r);
      phi += po->V(r);
    }
    
    outf << scientific 
         << setw(7) << rho << ' '
         << setw(12) << cbrt(4./rho) << ' '
         << setw(3) << layers.size() << ' '
         << setw(4) << np << ' '
         << setw(16) << rho_h << ' '
         << setw(16) << phi << ' '
         << setw(16) << po->F(rho_h) << ' '
         << setw(16) << po->F(rho_h) + 0.5*phi << endl;
  }
  
  outf.close();
  return;
} 

void write_subline_formnrg(Potential *po, sim_i *sim,
                           double Tmin, double Tmax, double step){
  double sub_density, zs, zsvi, dg = 0.;
  vector<vector<int*>> shells;
  int nm = 6;
  double *dgs = new double[nm];

  ofstream outf1,outf2;
  outf1.open("subline_db86.txt");
  outf2.open("formation_nrg_db86.txt");

  outf1 << setw(5) << "T" << ' '
        << setw(12) << "rho_sub" << endl;

  outf2 << setw(5) << "T" << ' '
        << setw(12) << "rho" << ' '
        << setw(12) << "Zs" << ' '
        << setw(12) << "ZsV1" << ' '
        << setw(12) << "ZsV2" << ' '
        << setw(12) << "ZsV3" << ' '
        << setw(12) << "ZsV4" << ' '
        << setw(12) << "ZsV5" << ' '
        << setw(12) << "ZsV6" << ' '
        << setw(12) << "dG1" << ' '
        << setw(12) << "dG2" << ' '
        << setw(12) << "dG3" << ' '
        << setw(12) << "dG4" << ' '
        << setw(12) << "dG5" << ' '
        << setw(12) << "dG6" << endl;

  for(double T = Tmin; T < Tmax; T += step){
    dg = 0.;
    sim -> T = T;
    sub_density = sublimation_density(po,sim);
    sim -> rho = sub_density;
    zs = config_integral_1p(po,sim);
    shells = get_shells(sim->from,sim->rho,nm);

    outf1 << setw(5) << T << ' '
          << setw(12) << setprecision(10) << sub_density << endl;

    outf2 << setw(5) << T << ' '
          << scientific << setw(12) << sub_density << ' '
          << setw(12) << zs << ' ';

    for(int i = 0; i < nm; i++){
      zsvi = config_integral_1p(po,sim,shells[i][0]);
      dg += -T*((double)(shells[i]).size())*log(zsvi/zs);
      outf2 << scientific << setw(12) << zsvi << ' ';
      dgs[i] = dg;
    }
    for(int i = 0; i < nm; i++){
      outf2 << scientific << setw(12) << dgs[i];
      if(i < nm-1) outf2 << ' ';
      else outf2 << endl;
    }
  }
  outf1.close();
  outf2.close();
  delete[] dgs;
  return;
}

int shared_ptcls(int *p1, int *p2, double rho, int nm){
  int res = 0;
  vector<vector<int*>> s1 = get_shells(p1,rho,nm);
  vector<vector<int*>> s2 = get_shells(p2,rho,nm);
  for(int i = 0; i < nm; i++){
    for(int j = 0; j < (int)(s1[i]).size(); j++){
      if(!(index_eq(s1[i][j],p2)) && (search_shells(s2,s1[i][j]) > 0)){
        print_indices(s1[i][j]);
        cout << " at shell " << i+1 << " from ";
        print_indices(p1);
        cout << " and shell " << search_shells(s2,s1[i][j])+1 << " from ";
        print_indices(p2);
        cout << endl;
        res+=1;
      }
    }
  }
  return res;
}

void hs_dg_2p(Potential *po, sim_i *sim){
  int vac_alpha[] = {0,1,0};
  int vac_beta[] = {-1,0,0};
  int vac_gamma[] = {-1,1,0};
  int vac_delta[] = {-1,1,1};
  int from2[] = {1,0,0};
  double zs = config_integral_1p(po,sim);
  double zsv1 = config_integral_1p(po,sim,from2);
  cout << "Zs: " << zs << endl << "ZsV1: " << zsv1 << endl;
  double alpha_correction = log(config_integral_2p(po,sim,from2,1,vac_alpha)/
                               (zsv1*zsv1));
  double beta_correction = log(config_integral_2p(po,sim,from2,1,vac_beta)/
                              (zsv1*zs));
  double gamma_correction = log(config_integral_2p(po,sim,from2,1,vac_gamma)/
                               (zsv1*zs));
  double delta_correction = log(config_integral_2p(po,sim,from2,1,vac_delta)/
                               (zsv1*zs));
  double ss_correction = log(config_integral_2p(po,sim,from2,1)/(zs*zs));
  double net_correction = -(sim->T)*(24*alpha_correction+12*beta_correction +
      48*gamma_correction+24*delta_correction-120*ss_correction);//+12*log(zsv1/zs));
  double p2p = pressure_2p(po,sim,6);
  cout << scientific << setprecision(10);
  //cout << "1p sublimation density is: " << sublimation_density(po,sim) << endl;
  //cout << "2p sublimation density is: " << sublimation_density_2p(po,sim,6) << endl;
  cout << "1p pressure at rho is " << pressure_1p(po,sim) << endl;
  cout << "2p pressure at rho is " << p2p << endl;
  cout << "1p dg is: " << formation_energy_1p(po,sim,6) << endl;
  cout << "alpha correction: " << alpha_correction << endl;
  cout << "beta  correction: " << beta_correction << endl;
  cout << "gamma correction: " << gamma_correction << endl;
  cout << "delta correction: " << delta_correction << endl;
  cout << "ss    correction: " << ss_correction << endl;
  cout << "net   correction: " << net_correction << endl;
  cout << "net       energy: " << p2p/(sim->rho)+net_correction << endl;
  return;
}

void hs_try2(Potential *po, sim_i *sim){
  double res = 0.0;
  double zs = config_integral_1p(po,sim);
  vector<vector<int*>> nn = get_shells(sim->from,sim->rho,1);
  vector<vector<int*>> max_shells = get_shells(sim->from,sim->rho,20);
  for(int i = 0; i < (int)nn[0].size(); i++){
    int *from2 = nn[0][i];
    double zsv1 = config_integral_1p(po,sim,from2);
    vector<vector<int*>> nnn = get_shells(nn[0][i],sim->rho,1);
    for(int j = 0; j < (int)nnn[0].size(); j++){
      if(!index_eq(sim->from,nnn[0][j])){
        if(search_shells(max_shells,nnn[0][j]) == 0)
          res+=log(config_integral_2p(po,sim,from2,1,nnn[0][j])/(2*zsv1*zsv1));
        else
          res += log(config_integral_2p(po,sim,from2,1,nnn[0][j])/(zsv1*zs));
      }
    }
  }
  double dg1p = formation_energy_1p(po,sim,6);
  cout << scientific << setprecision(10);
  cout << "1p pressure at rho is: " << pressure_1p(po,sim) << endl;
  cout << "2p pressure at rho is: " << pressure_2p(po,sim,6) << endl;
  cout << "1p dg is: " << dg1p << endl;
  cout << "net correction is: " << -(sim->T)*res << endl;
  cout << "net result 6 shell is: " << -(sim->T)*(res+dg1p) << endl;;
  return;
}

int main(int argc, char* argv[]){
  if(argc != 3){
    cout << "USAGE: ./run T rho" << endl;
  }

  int from[] = {0,0,0};
  sim_i *sim = new sim_i;
  sim->from = from;
  sim->T = atof(argv[1]);
  sim->rho = atof(argv[2]);
  sim->asd = ASD(sim->rho);
  
  LJ potential(2*cbrt(sqrt(2)));///6.0);
  //eamf85 potential(ADATA,EDATA_F);
  //eamdb86 potential(ADATA,EDATA_DB);
  //eamm99 potential(ADATA,EDATA_M);
  
  //free_energy_2p(&potential,sim,1);
  double subd = sublimation_density(&potential,sim);
  sim -> rho = subd;
  sim -> asd = ASD(sim->rho);
  //hs_try2(&potential,sim);
  hs_dg_2p(&potential,sim);
  //cout << vac_nrgvd(&potential, sim, 1)[0] << endl;
  //cout << setprecision(10) << sublimation_density_2p(&potential,sim,6) << endl;
  
  //for(double T = 0.1; T < 1.0; T += 0.1)
  //  cout << T << ' ' << vac_nrgvd(&potential,(sim->T=T,sim->from=from,sim),1)[0] << endl;

  //for(int nm = 8; nm < 11; nm++)
  //  cout << nm << ' ' << vac_nrgvd(&potential,(sim->from=from,sim),nm)[0] << endl;

  return 0;
}
