import subprocess
import numpy as np
import matplotlib.pyplot as plt

OUTNAME = "../output_files/eos_analyze_out.txt"
OUTNAME2P = "../output_files/f_v_rho_2ptcl_out.txt"

def plot_output():
  rho, zs, vdh = np.loadtxt(OUTNAME,unpack = True)
  T2, rho2, f1, c1, c2, c3, c4, c5, c6, f2 = np.loadtxt(OUTNAME2P, unpack=True)
  p1, = plt.plot(rho,zs,marker = 'o')
  p2, = plt.plot(rho,vdh,marker = 'o')
  p3, = plt.plot(rho2,f2,marker = 'o')
  plt.legend([p1,p2,p3],["1st Order", "Van der Hoef", "2nd Order"],loc=4)
  plt.xlabel("$\\rho$")
  plt.ylabel("F")
  plt.title("Helmholz Free Energy vs. Density")
  plt.show()

def main():
  outf = open("../output_files/eos_analyze_out.txt",'w')
  T = 1.0
  rhos = np.arange(0.2,2.0,.1)

  for rho in rhos:
    cmd1 = ['../src/get_energy', str(T), str(rho)]
    cmd2 = ['../test_scripts/eos_tests/run', str(T), str(rho)]
    
    p1 = subprocess.Popen(cmd1,stdout = subprocess.PIPE,\
                          stderr=subprocess.PIPE,stdin=subprocess.PIPE)
    out1, err1 = p1.communicate()

    p2 = subprocess.Popen(cmd2,stdout = subprocess.PIPE,\
                          stderr=subprocess.PIPE,stdin=subprocess.PIPE)
    out2, err2 = p2.communicate()

    outf.write(str(rho)+'\t'+out1[:len(out1)-1]+'\t'+out2)

  outf.close()
  plot_output()

plot_output()
