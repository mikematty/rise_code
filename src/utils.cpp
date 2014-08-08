#include "utils.hpp"

bool index_eq(int *p1, int *p2){
  return ((p1[0] == p2[0]) && (p1[1] == p2[1]) && (p1[2] == p2[2]));
}

void print_indices(int *f){
  cout << "(" << f[0] << "," << f[1] << "," << f[2] << ")";
  return;
}


double fcc_raw(double *c, int *from, int *to, double asd, double *c2){
  
  //want to treat (l1,m1,n1) as (0,0,0)
  double l = (double)(to[0]-from[0]);
  int m = to[1]-from[1];
  int n = to[2]-from[2];
  
  double nmod = (double)((n%3+3)%3);
  double mmod = (double)(((to[1]%2+2)%2)-((from[1]%2+2)%2));
  double factor = nmod*(-1.5*nmod+2.5)*0.5;
  if(((to[2]%3+3)%3)*((from[2]%3+3)%3) > 0) factor *= -2;
  double s3 = sqrt(3);
  double x = (factor+l+mmod*0.5)*asd-c[0];
  double y = (-factor/s3+m*s3/2.)*asd-c[1];
  double z = (n*sqrt(2.)/s3)*asd-c[2];

  //if the second particle is off center need to add in its coordinates
  if(c2 != NULL){
    x += c2[0];
    y += c2[1];
    z += c2[2];
  }

  return (x*x+y*y+z*z);
}

double fcc_distance(double *c, int *from, int *to, double asd, double *c2){
  //want to treat (l1,m1,n1) as (0,0,0)
  double l = (double)(to[0]-from[0]);
  int m = to[1]-from[1];
  int n = to[2]-from[2];
  
  double nmod = (double)((n%3+3)%3);
  double mmod = (double)(((to[1]%2+2)%2)-((from[1]%2+2)%2));
  double factor = nmod*(-1.5*nmod+2.5)*0.5;
  if(((to[2]%3+3)%3)*((from[2]%3+3)%3) > 0) factor *= -2;
  double s3 = sqrt(3);
  double x = (factor+l+mmod*0.5)*asd-c[0];
  double y = (-factor/s3+m*s3/2.)*asd-c[1];
  double z = (n*sqrt(2.)/s3)*asd-c[2];

  //if the second particle is off center need to add in its coordinates
  if(c2 != NULL){
    x += c2[0];
    y += c2[1];
    z += c2[2];
  }

  return sqrt(x*x+y*y+z*z);
}

bool query_ws(int *f, double rho, double *pt){
  bool res = true;
  //vector<vector<int*>> nn = get_shells(f, rho, 1); // nearest neighbors
  int nn[12][3] = {{0,-1,0},{-1,1,0},{1,-1,-1},{1,-1,0},{-1,1,1},{-1,0,1},
                   {-1,0,0},{0,1,0},{0,0,-1},{0,0,1},{1,0,-1},{1,0,0}};
  double mag = sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]); //distance to pt
  for(int i = 0; i < 12; i++){//(int)(nn[0]).size(); i++){
    //by construction of ws cell, if pt closer to any neighbor, not in ws cell
    if (fcc_distance(pt, f, nn[i], ASD(rho)) < mag) res = false;
    //delete[] nn[0][i];
  }
  return res;
}

int search_shells(vector<vector<int*>> shells, int *target){
  for(int i = 0; i < (int)shells.size(); i++){
    for(int j = 0; j < (int)(shells[i]).size(); j++){
      if (index_eq(target,shells[i][j])) return i; //found it!
    }
  }
  return -1; //couldn't find target :(
}

vector<int*> find_contributors(int *f, double rho, double max_dist){
  vector<int*> retv(0);
  int* test_indices;
  int x_max, y_max, z_max;
  double *coords = new double[3](); //want to check from center
  int *t = new int[3]; //index vector for sphere to check (t is for to)

  //build rectangle of ptcls where farthest on each side is w/in max_dist
  for(int i = 0; i < 3; i++){
    t[i] = f[i]+1;
    t[(i+1)%3] = f[(i+1)%3];
    t[(i+2)%3] = f[(i+2)%3];
    while(fcc_distance(coords,f,t,ASD(rho)) < max_dist+0.5*ASD(rho)) t[i] += 1;
    if(i == 0) x_max = t[i];
    else if(i == 1) y_max = t[i];
    else z_max = t[i];
  }

  //add those in the rectangle to the result that are actually within max_dist
  for(int i = f[0]-(x_max-f[0]); i <= x_max; i++){
    for(int j = f[1]-(y_max-f[1]); j <= y_max; j++){
      for(int k = f[2]-(z_max-f[2]); k <= z_max; k++){
        test_indices = new int[3];
        test_indices[0] = i;
        test_indices[1] = j;
        test_indices[2] = k;
        //push the allowed indices onto the resultant set of contributors
        if((fcc_distance(coords,f,test_indices,ASD(rho)) < 
            max_dist+0.5*ASD(rho)) && !index_eq(test_indices,f))
          retv.push_back(test_indices);
      } //i loop
    } //j loop
  } //k loop
  //clean up
  delete [] coords;
  delete [] t;
  return retv;
}

vector<vector<int*>> get_shells(int *f, double rho, double max_dist){
  vector<int*> cs = find_contributors(f, rho, max_dist);
  double *origin = new double[3]();
  //lambda function to order indices based on distance from f
  auto ordering_func = [=](int* i, int* j)->bool {
    double idist = fcc_distance(origin,f,i,ASD(rho));
    double jdist = fcc_distance(origin,f,j,ASD(rho));
    return (idist < jdist);
  };
  sort(cs.begin(),cs.end(),ordering_func); //sort all relevant particles

  //next separate into shells
  double d; //distance to particle
  double prev_d = 0.; //distance to previously visited particle
  int layer = 0; //shell index
  //vector to hold particles within current layer
  vector<int*>* cur_layer = new vector<int*>;
  vector<vector<int*>> shells;
  for(int i = 0; i < ((int)cs.size()); i++){
    d = fcc_distance(origin,f,cs[i],ASD(rho));  
    //check if current distance is not equal to previous distance
    //use finite difference instead of equivalence b/c floating point error
    if(fabs(d-prev_d) > 1e-3){
      if(layer > 0){
        shells.push_back(*cur_layer); //add shell to result
        cur_layer = new vector<int*>;
      }
      layer += 1; //now we're on the next shell
    }
    cur_layer->push_back(cs[i]); //add particle to shell
    prev_d = d;
  }
  shells.push_back(*cur_layer);
  delete[] origin; //clean up
  return shells;
}

vector<vector<int*>> get_shells(int *f, double rho, int nm){
  int overkill_index[] = {f[0]+nm,0,0};
  double *origin = new double[3]();
  double max_dist = fcc_distance(origin,f,overkill_index,ASD(rho));
  //cout << "max_dist is " << max_dist << endl;
  vector<int*> cs = find_contributors(f, rho, max_dist);
  auto ordering_func = [=](int* i, int* j)->bool {
    double idist = fcc_distance(origin,f,i,ASD(rho));
    double jdist = fcc_distance(origin,f,j,ASD(rho));
    return (idist < jdist);
  };
  sort(cs.begin(),cs.end(),ordering_func);

  //next separate into shells
  double d;
  double prev_d = 0.;
  int layer = 0;
  vector<int*>* cur_layer = new vector<int*>;
  vector<vector<int*>> shells;
  for(int i = 0; i < ((int)cs.size()); i++){
    d = fcc_distance(origin,f,cs[i],ASD(rho));  
    //print_indices(cs[i]); cout << " is " << d << " from "; print_indices(f); cout << endl;
    if(fabs(d-prev_d) > 1e-3){
      if(layer > 0){
        shells.push_back(*cur_layer);
        cur_layer = new vector<int*>;
      }
      layer += 1;
      if(layer > nm) break;
    }
    cur_layer->push_back(cs[i]);
    prev_d = d;
  }
  //if(nm == 1) shells.push_back(*cur_layer);
  delete[] origin;
  return shells;
}

vector<int*> shell_intersect(vector<vector<int*>> shells1, 
                             vector<vector<int*>> shells2) {
  vector<int*> retv(0);
  for(int i = 0; i < (int)shells1.size(); i++){
    for(int j = 0; j < (int)(shells1[i]).size(); j++){
      if(search_shells(shells2,shells1[i][j]) >= 0) retv.push_back(shells1[i][j]);
    }
  }
  return retv;
}

void print_ptcl_table(int *f, double rho, vector<vector<int*>> shells){
  int ptcl = 1;
  double origin[] = {0.,0.,0.};

  ofstream outf;
  outf.open("../output_files/ptcl_table.txt");
  outf << setw(6) << "ptcl" << '\t'
       << setw(5) << "layer" << '\t'
       << setw(4) << "pnl" << '\t'
       << setw(8) << "d" << '\t'
       << "(" << setw(3) << "l"
       << "," << setw(3) << "m"
       << "," << setw(3) << "n" << ")" << '\t'
       << "(" << setw(10) << "x"
       << "," << setw(10) << "y"
       << "," << setw(10) << "z" << ")"
       << endl;

  for(int i = 0; i < (int)shells.size(); i++){
    for(int j = 0; j < (int)(shells[i]).size(); j++){
      outf << setw(6) << ptcl << '\t'
           << setw(5) << i+1 << '\t'
           << setw(4) << j+1 << '\t'
           << setw(8) << fcc_distance(origin,f,shells[i][j],ASD(rho)) << '\t'
           << "(" << setw(3) << shells[i][j][0] << ","
           << setw(3) << shells[i][j][1] << ","
           << setw(3) << shells[i][j][2] << ")" << '\t'
           << endl;
      ptcl += 1;
    }
    outf << endl;
  }
  outf.close();
  return;
}

