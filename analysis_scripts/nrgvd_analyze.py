import sys, numpy as np, matplotlib.pyplot as plt

def main():
  if(len(sys.argv) != 4):
    print "USAGE: python nrgvd_analyze.py T rho nm"
    return
  T = sys.argv[1]
  rho = sys.argv[2]
  nm = sys.argv[3]

  fname = "../output_files/vac_int_data/vac_nrgvd_nm"+nm+\
          "_t"+T.replace('.','-')+'_rho'+rho.replace('.','-')+'.txt'

  inf = open(fname)
  data = inf.readlines()

  seps = [np.sqrt(float(x.split()[2])/2.)*(4./float(rho))**(1./3) for x in data]
  dfs = [float(x.split()[5]) for x in data]

  plt.plot(seps,dfs,marker='o',linestyle='')
  plt.xlabel("Separation Distance")
  plt.ylabel("$\Delta F$")
  plt.title("T = "+T+", rho = "+rho)

  plt.show()

if __name__ == "__main__":main()
