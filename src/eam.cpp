#include "eam.hpp"

/******************** EAM BASE CLASS IMPLEMENTATION **************************/
eam::eam(){}; //empty default constructor for a.b.c.

//constructor to be called from derived classes
eam::eam(string atom_data){
  read_atom_data(atom_data);
  get_R_coeffs();
};

eam::~eam(){};

double eam::exptential_1p(vector<int*> cs, sim_i *sim, double *coords){
  double rho_h = 0.; //net electronic density
  double phi = 0.;  //net pair potential
  double r;
  for(int i = 0; i < (int)cs.size(); i++){
    r = fcc_distance(coords,sim->from,cs[i],sim->asd);
    rho_h += rho_a(r);
    phi += V(r);
  }
  //return the potential
  return exp(-BETA(sim->T)*(F(rho_h) + 0.5*phi));
}

double eam::exptential_2p(vector<int*> cs1, vector<int*> cs2,
                          sim_i *sim, int *from2, double *coords){
  //net pair potentials and electronic denisties for both particles
  double rho_h1 = 0., rho_h2 = 0.;
  double phi1 = 0., phi2 = 0.;
  double r;
  //split into 1 particle (x,y,z) coords for each particle
  double coords1[] = {coords[0], coords[1], coords[2]};
  double coords2[] = {coords[3], coords[4], coords[5]};
  //both sets of contributors should be same size so only need 1 loop
  for(int i = 0; i < (int)cs1.size(); i++){
    //if at 2nd particle, need to account for its coordinates
    if(index_eq(cs1[i],from2)) 
      r = fcc_distance(coords1,sim->from,cs1[i],sim->asd,coords2);
    else r = fcc_distance(coords1,sim->from,cs1[i],sim->asd);

    rho_h1 += rho_a(r);
    phi1 += V(r);
    //if at 1st particle, need to account for its coordinates
    if(index_eq(cs2[i],sim->from))
      r = fcc_distance(coords2,from2,cs2[i],sim->asd,coords1);
    else r = fcc_distance(coords2,from2,cs2[i],sim->asd);

    rho_h2 += rho_a(r);
    phi2 += V(r);
  }
  return exp(-BETA(sim->T)*((F(rho_h1)+0.5*phi1)+
                            (F(rho_h2)+0.5*phi2)));
}

double eam::rhoa_s(double r){
  double res = 0.;
  //just compute based on form in paper
  for(int i = 0; i < (int)zeta_s.size(); i++)
    res += Cs[i]*Rs[i]*pow(r,i/2)*exp(-zeta_s[i]*r);
  return (res*res)/(4*M_PI);
}

double eam::rhoa_d(double r){
  double res = 0.;
  for(int i = 0; i < (int)zeta_d.size(); i++)
    res += Cd[i]*Rd[i]*pow(r,2+i/2)*exp(-zeta_d[i]*r);
  return (res*res)/(4*M_PI);
}

void eam::get_R_coeffs(){
  int n_i,fact;
  //precompute these for efficiency reasons
  //just based on form in paper
  for(int i = 0; i < (int)zeta_s.size(); i++){
    n_i = i/2 + 1;
    fact = 1;
    for(int j = 1; j <= 2*n_i; j++) fact *= j; 
    Rs.push_back(pow(2*zeta_s[i],n_i+0.5)/sqrt((double)fact));
  }
  for(int i = 0; i < (int)zeta_d.size(); i++){
    n_i = i/2+3;
    fact = 1;
    for(int j = 1; j <= 2*n_i; j++) fact *= j;
    Rd.push_back(pow(2*zeta_d[i],n_i+0.5)/sqrt((double)fact));
  }
  return;
}

void eam::read_atom_data(string in_name){
  //just read input data in specific format
  //example of format given in rise_code/input_files/Ni_atom.dat
  char ch = ' ';
  int i = 1;
  double zeta, c;

  fstream inf;
  inf.open(in_name, ios::in);
  
  while(ch != ':') inf >> ch;
  while(i != 0) {
    inf >> i >> zeta >> c;
    zeta_s.push_back(zeta);
    Cs.push_back(c);
  }

  ch = ' ';
  i = 1;
  while(ch != ':') inf >> ch;
  while(i != 0) {
    inf >> i >> zeta >> c;
    zeta_d.push_back(zeta);
    Cd.push_back(c);
  }

  ch = ' ';
  while(ch != ':') inf >> ch;
  inf >> rho0 >> k0 >> kp0;

  inf.close();
  return;
}

/******************* END EAM BASE CLASS IMPLEMENTATION ***********************/

/****************** DAW/BASKES 86 EAM POTENTIAL IMPLEMENTATION ***************/
eamdb86::eamdb86(string atom_data, string eam_data): eam(atom_data) {
  read_eam_data(eam_data);
  get_knots();
  Fpp_knots = spline_d2(rho_a_knots,F_knots,DMAX,DMAX);
  max_dist = 0.1;
  while(fabs(F(rho_a(max_dist)) + V(max_dist)) > 1e-4) max_dist += 0.1;
}

eamdb86::~eamdb86(){};

double eamdb86::rho_a(double r){
  return ns*rhoa_s(r) + (n-ns)*rhoa_d(r);
}

double eamdb86::V(double r) {
  double r_nu = (nu == 2) ? r*r : r; //compute r^nu
  double Z = z0 * (1 + beta*r_nu) * exp(-alpha*r); //equationf rom db86
  return Z*Z*FtoDB /r;
}

double eamdb86::F(double rho) {
  int n = (int)rho_a_knots.size();
  if(rho < 0.) return 0.;
  else if(rho > rho_a_knots[n-1]){
    //extrapolate a value for F in this case (linear extrapolation)
    return F_knots[n-1] + ((F_knots[n-1]-F_knots[n-2])/
           (rho_a_knots[n-1]-rho_a_knots[n-2]) + 
           (rho_a_knots[n-1]-rho_a_knots[n-2]) * Fpp_knots[n-2] / 6.) *
           (rho - rho_a_knots[n-1]);
  }
  //otherwise use spline interpolation to find a value for F
  else return splint(rho_a_knots,F_knots,Fpp_knots,rho);
}

void eamdb86::read_eam_data(string in_name){
  //example of data format in
  //rise_code/input_files/Ni_EAM_DB86.dat
  char ch = ' ';
  fstream inf;
  inf.open(in_name, ios::in);

  while(ch != ':') inf >> ch;
  inf >> n >> ns;

  ch = ' ';
  while(ch != ':') inf >> ch;
  inf >> z0 >> alpha >> beta >> nu;

  ch = ' ';
  while(ch != ':') inf >> ch;
  inf >> a0 >> Esub0 >> B;

  inf.close();
  return;
}

void eamdb86::get_knots(){
  double a_dev = 10.; //lattice constant deviation
  double a_res = 20.; //lattice constant range resolution
  double a_max = (1. + a_dev / 10.) * a0; //maximum lattice constant
  double a_min = (1. - a_dev / 100.) * a0; //minimum lattice constant
  double rc = 8.; //temporary cutoff range
  int f[] = {0,0,0};
  double origin[] = {0.,0.,0.};
  double omega = 4./(a0*a0*a0); //equilibrium atomic volume
  double rho,r,phi,rho_bar;
  double a_star, Esub;
  vector<vector<int*>> shells;
  //know at 0 eletronic density, embedding function should be 0
  rho_a_knots.push_back(0.);
  F_knots.push_back(0.);
  for(double a = a_max; a >= a_min; a -= (a_max-a_min)/a_res){
    rho = 4./(a*a*a); //particle density
    shells = get_shells(f, rho, rc);
    phi = rho_bar = 0.;
    //sum up individual contributions to e- density and pair potential
    for(int i = 0; i < (int)shells.size(); i++){
      r = fcc_distance(origin,f,shells[i][0],ASD(rho));
      phi += (double)(shells[i]).size()*V(r);
      rho_bar += (double)(shells[i]).size()*rho_a(r);
    }
    a_star = (a/a0-1.)/sqrt(Esub0/(9. * B * omega));
    Esub = -Esub0*(1+a_star)*exp(-a_star); //sublimation energy
    rho_a_knots.push_back(rho_bar);
    F_knots.push_back(Esub-0.5*phi); //equation from db86
  }
  return;
}
    
/*********** END DAW/BASKES 86 EAM POTENTIAL IMPLEMENTATION ******************/

/*************** FOILES 85 EAM POTENTIAL IMPLEMENTATION **********************/

eamf85::eamf85(string atom_data, string eam_data): eam(atom_data) {
  read_eam_data(eam_data);
  Fpp_knots = spline_d2(rho_knots,F_knots,DMAX,DMAX);
  max_dist = 0.1;
  while(fabs(F(rho_a(max_dist))+V(max_dist)) > 1e-4) max_dist += 0.1;
}

eamf85::~eamf85(){};

double eamf85::rho_a(double r){
  return ns*rhoa_s(r) + (n-ns)*rhoa_d(r);
}

double eamf85::V(double r){
  if(r > ZRC) return 0.;
  return (Za1*pow(ZRC-r,3)+Za2*pow(ZRC-r,4))*FtoDB/r;
}

double eamf85::F(double rho) {
  int n = (int)rho_knots.size();
  if(rho < 0.) return 0.;
  else if(rho > rho_knots[n-1]){
   //linearly extrapolate a value for F in this case
    double un = (F_knots[n-1]-F_knots[n-2])/(rho_knots[n-1]-rho_knots[n-2])+
                (rho_knots[n-1]-rho_knots[n-2])*Fpp_knots[n-2]/6.;
    return F_knots[n-1]+un*(rho-rho_knots[n-1]);
  }
  //otherwise use spline interpolation to find a value for F
  else return splint(rho_knots,F_knots,Fpp_knots,rho);
}

void eamf85::read_eam_data(string in_name){
  //for an example of input data format, see
  //rise_code/input_files/Ni_EAM_F85.dat
  char ch = ' ';
  double a,b;
  int i = 1;
  fstream fin;
  fin.open(in_name, ios::in);

  while(ch != ':') fin >> ch;
  fin >> n >> ns;

  ch = ' ';
  while(ch != ':') fin >> ch;
  while(i != 0) {
    fin >> i >> a >> b;
    rho_knots.push_back(a);
    F_knots.push_back(b);
  }

  ch = ' ';
  while(ch != ':') fin >> ch;
  fin >> Za1 >> Za2 >> ZRC;

  fin.close();
  return;
}

/************** END FOILES 85 EAM POTENTIAL IMPLEMENTATION *******************/

/*************** MISHIN 99 EAM POTENTIAL IMPLEMENTATION **********************/
eamm99::eamm99(string atom_data, string eam_data): eam(atom_data) {
  read_eam_data(eam_data);
  V_cutoff = r_knots[r_knots.size()-1]; //after this point V is 0
  rhoa_cutoff = r_knots[r_knots.size()-1-n_empty]; //rhoa is 0 after here

  vector<double> r_knots_temp = r_knots;
  r_knots_temp.erase(r_knots_temp.begin()+r_knots_temp.size()-n_empty,
                     r_knots_temp.begin()+r_knots_temp.size());
  rhoapp_knots = spline_d2(r_knots_temp,rhoa_knots,DMAX,DMAX);
  for(int i = 0; i < n_empty; i++) rhoapp_knots.push_back(0.);

  Vpp_knots = spline_d2(r_knots,V_knots,DMAX,DMAX);
  Fpp_knots = spline_d2(rho_knots,F_knots,DMAX,DMAX);
}

eamm99::~eamm99(){};

double eamm99::rho_a(double r){
  //if out negative or too large, just return 0
  if((r < 0) || (r > r_knots[r_knots.size()-1])) return 0.;
  else if(r < r_knots[0]){
    //linearly extrapolate a value for rho_a
    double res = rhoa_knots[0]-(r_knots[0]-r)*
                 ((rhoa_knots[1]-rhoa_knots[0])/(r_knots[1]-r_knots[0]) - 
                  (r_knots[1]-r_knots[0])*rhoapp_knots[1]/6.);
    return ((res > 0.) ? res : 0.);
  }
  //if in range, use cubic spline interpolation
  else return splint(r_knots,rhoa_knots,rhoapp_knots,r);
}

double eamm99::V(double r){
  if((r < 0) || (r > r_knots[r_knots.size()-1])) return 0.;
  else if(r < r_knots[0]){
    //etrapolate a value for V
    double un = (V_knots[1] - V_knots[0])/(r_knots[1]-r_knots[0]) - 
                (r_knots[1]-r_knots[0])*Vpp_knots[1]/6.0;
    return V_knots[0]-un*(r_knots[0]-r);
  }
  //if in range, use cubic spline interpolation
  else return splint(r_knots,V_knots,Vpp_knots,r);
}

double eamm99::F(double rho){
  int n = (int)rho_knots.size();
  if(rho < 0.) return 0.; //sure
  else if(rho >= rho_knots[n-1]){
   //extrapolate a value for F in this case
    return F_knots[n-1] + ((F_knots[n-1]-F_knots[n-2])/
           (rho_knots[n-1]-rho_knots[n-2]) + (rho_knots[n-1]-rho_knots[n-2]) * 
           Fpp_knots[n-2] / 6.) * (rho - rho_knots[n-1]);
  }
  //otherwise use spline interpolation to find a value for F
  else return splint(rho_knots,F_knots,Fpp_knots,rho);
}

void eamm99::read_eam_data(string in_name){
  //for an example of data format, see
  //rise_code/input_files/Ni_EAM_M99.dat
  double a,b,c,d,e;
  char ch = ' ';
  int i = 1;

  fstream fin;
  fin.open(in_name, ios::in);

  while(ch != ':') fin >> ch;
  while(i != 0){
    fin >> i >> a >> b >> c >> d >> e;
    r_knots.push_back(a);
    V_knots.push_back(b);
    rhoa_knots.push_back(c);
    rho_knots.push_back(d);
    F_knots.push_back(e);
  }

  ch = ' ';
  while(ch != ':') fin >> ch;
  fin >> n_empty;

  fin.close();
  max_dist = r_knots[r_knots.size()-1];
  return;
}

/************* END MISHIN 99 EAM POTENTIAL IMPLEMENTATION ********************/
