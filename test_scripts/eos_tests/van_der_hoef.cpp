#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#define T_tr 0.687

using namespace std;

namespace hoef{
  const double a[3][4] = {{-8.2151768,12.070686,-6.6594615,1.3211582},
                          {13.404069,-20.632066,11.564825,-2.3064801},
                          {-5.5481261,8.8465978,-5.0258631,1.0070066}};

  const double b_n[4] = {69.833875,-132.86963,97.438593,-25.848057};

  const double c_n[5] = {0.,0.,-14.45392093,0.,6.065940096};

  double beta_a_ex(double temperature, double density);
}

double hoef::beta_a_ex(double T, double rho){
  double C, u_stat, U_ah = 0., b = 0.;
  double beta = 1./T;

  C = -19.4503982 - 8.89103855 * T_tr + 4.68985418 * T_tr * T_tr;
  u_stat = c_n[2] * rho * rho + c_n[4] * pow(rho, 4);


  for(int n = 0; n <= 2; n++){
    for(int m = 0; m <= 2; m++){
      U_ah -= a[n][m] * pow(rho, n) * pow(beta, -m-1) / ((double)(m+1));
    }
  }

  for(int n = 0; n <= 3; n++){
    b += b_n[n] * pow(rho, n+1) / ((double)(n+1));
  }

  return (C+beta*u_stat+1.5*log(beta)+U_ah+b);
}

int main(int argc, char*argv[]){
  double T = atof(argv[1]);
  double rho = atof(argv[2]);
  cout << setprecision(10) << hoef::beta_a_ex(T,rho)+rho*(log(rho)-1) << endl;
  return 0;
}
