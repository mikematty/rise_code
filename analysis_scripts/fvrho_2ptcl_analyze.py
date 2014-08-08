from subprocess import *
import numpy as np
import matplotlib.pyplot as plt

OUTNAME = "../output_files/f_v_rho_2ptcl_ljshort_out.txt"

def main():
  outf = open(OUTNAME,'w')
  outf.write("#T\trho\t")
  outf.write("first order\tcorrection 1\tcorrection 2\tcorrection 3\t")
  outf.write("correction4\tcorrection 5\tcorrection 6\tsecond order \n")
  
  #LJ
  Ts = [0.3, 0.6, 0.9, 1.2]
  rhos = np.arange(0.95,1.17,.02)
  
  #EAM
  #Ts = [.09,.11,.13,.15,.17]
  #rhos = np.linspace(.08,.1,11)

  for T in Ts:
    for rho in rhos:
      cmd = ['../src_v2/run',str(T),str(rho)]
      p = Popen(cmd,stdout=PIPE,stderr=PIPE,stdin=PIPE)
      out,err = p.communicate()
      terms = out.split()
      terms.sort()
      terms = map(lambda x: x[2:],terms)
      
      outf.write(str(T)+'\t'+str(rho)+'\t')
      print T,rho,
      for n in terms:
        #print n, '\t',
        outf.write(n+'\t')
      outf.write('\n')
      print

  outf.close()

if __name__ == "__main__":main()
