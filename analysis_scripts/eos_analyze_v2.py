from subprocess import PIPE,Popen
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

PO = "ljshort"
OUTNAME = "../output_files/eos_analyze_ljshort_out.txt"
INNAME = "../output_files/f_v_rho_2ptcl_ljshort_out.txt"

def plot_fe_output(n):
  Ts,rhos,f1,f2,vdh = np.loadtxt(OUTNAME,unpack=True)
  sub_plot = 221
  for i in xrange(0,len(rhos),n):
    plt.subplot(sub_plot)
    s = slice(i,i+n)
    p1, = plt.plot(rhos[s],f1[s],marker='o',linewidth=3)
    p2, = plt.plot(rhos[s],f2[s],marker='>',linewidth=3,linestyle='-.')
    p3, = plt.plot(rhos[s],Ts[i]*vdh[s],marker='d',linewidth=3,linestyle=':')
    plt.xlabel('$\\rho$',fontsize=20)
    plt.ylabel('F',fontsize=15)
    plt.title('F vs $\\rho$ for T = '+str(Ts[i]))
    plt.legend([p1,p2,p3],['$1^{st}$ order','$2^{nd}$ order','Van der Hoef'],\
               loc=2)
    sub_plot += 1
  plt.show()
  return
  
def cubic(x,*args):
  a,b,c,d = args
  return a*x**3+b*x**2+c*x+d

def plot_p_output(n):
  Ts,rhos,f1,f2,vdh = np.loadtxt(OUTNAME,unpack=True)
  plot_mrho = []
  plot_mhoef = []
  plot_mcalc = []
  nmostafa = 400
  sub_plot = 221
  #fn = lambda xs:(lambda a,b,c:[a*xs[i]**2+b*xs[i]+c for i in xrange(len(xs))])
  #quad = fn(rhos)
  quad = lambda a,b,c:[(rhos[i]**2)*(3*a*rhos[i]**2+2*b*rhos[i]+c)\
                       for i in xrange(len(rhos))]

  for i in xrange(0,len(rhos),n):
    s = slice(i,i+n)
    mslice = slice((i/n)*nmostafa,((i/n)+1)*nmostafa)
    plt.subplot(sub_plot)
    popt1,_ = curve_fit(cubic,rhos[s],f1[s],p0=[1,1,1,1])
    popt2,_ = curve_fit(cubic,rhos[s],f2[s],p0=[1,1,1,1])
    poptvdh,_ = curve_fit(cubic,rhos[s],Ts[i]*vdh[s],p0=[1,1,1,1])
    p1, = plt.plot(rhos[s],quad(popt1[0],popt1[1],popt1[2])[s],\
                   marker = 'o', linewidth = 3)
    p2, = plt.plot(rhos[s],quad(popt2[0],popt2[1],popt2[2])[s],\
                   marker = '>', linewidth = 3, linestyle = '-.')
    p3, = plt.plot(rhos[s],quad(poptvdh[0],poptvdh[1],poptvdh[2])[s],\
                   marker = 'd', linewidth = 3, linestyle = ':')
    plt.xlabel('$\\rho$',fontsize = 20)
    plt.ylabel('P',fontsize = 15)
    plt.title('P vs $\\rho$ for T = '+str(Ts[i]))
    plt.legend([p1,p2,p3],['$1^{st}$ order','$2^{nd}$ order','Van der Hoef'],\
               loc = 2)
    sub_plot += 1
  plt.show()
  return

def plot_subline(n):
  Ts,rhos,f1,f2,fhoef = np.loadtxt(OUTNAME,unpack=True)
  plot_mrho = []
  plot_mhoef = []
  plot_mcalc = []
  subd1 = []
  subd2 = []
  tpts = []
  quad = lambda a,b,c:[(rhos[i]**2)*(3*a*rhos[i]**2+2*b*rhos[i]+c)\
                       for i in xrange(len(rhos))]
  quad_psolve = lambda a,b,c: (-b+np.sqrt(b*b-4*a*c))/(2*a)
  
  for i in xrange(0,len(rhos),n):
    s = slice(i,i+n)
    popt1,_ = curve_fit(cubic,rhos[s],f1[s],p0=[1,1,1,1])
    popt2,_ = curve_fit(cubic,rhos[s],f2[s],p0=[1,1,1,1])
    subd1 += [quad_psolve(popt1[0]*rhos[i]**2*3,popt1[1]*rhos[i]**2*2,\
                          popt1[2]*rhos[i]**2)]
    subd2 += [quad_psolve(popt2[0]*rhos[i]**2*3,popt2[1]*rhos[i]**2*2,\
                          popt2[2]*rhos[i]**2)]
    tpts += [Ts[i]]
    print tpts[-1],'\t',subd1[-1],'\t',subd2[-1]

  p1, = plt.plot(tpts,subd1,linewidth=3,marker='o')
  p2, = plt.plot(tpts,subd2,linewidth=3,linestyle='-.',marker='>')
  plt.xlabel('$T(K)$',fontsize=20)
  plt.ylabel('$\\rho(g/cm^3)$',fontsize=20)
  plt.title("Sublimation Density vs. Temperature")
  plt.legend([p1,p2],['$1^{st}$ order','$2^{nd}$ order'],loc=1)
  plt.show()
  return

def main():
  outf = open(OUTNAME,'w')
  Ts = [0.3, 0.6, 0.9, 1.2]
  rhos = np.arange(0.95,1.17,0.02)
  vdh_pts = []

  for T in Ts:
    for rho in rhos:
      cmd = ['../test_scripts/eos_tests/run',str(T),str(rho)]
      p = Popen(cmd,stdout=PIPE,stderr=PIPE,stdin=PIPE)
      out,err = p.communicate()
      vdh_pts += [out]

  T_in,rho_in,f1,f2 = np.loadtxt(INNAME,usecols=[0,1,2,9],unpack = True)
  
  for i in xrange(len(T_in)):
    outf.write(str(T_in[i])+'\t'+str(rho_in[i])+'\t'+str(f1[i])+'\t'+\
               str(f2[i])+'\t'+vdh_pts[i])
  outf.close()
  plot_fe_output(len(rhos))
  plot_p_output(len(rhos))
  plot_subline(len(rhos))

main()
