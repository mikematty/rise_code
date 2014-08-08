import numpy as np
import matplotlib.pyplot as plt

def main():
  ts,ps,rhos = np.loadtxt("../output_files/p_v_rho.txt", unpack = True)
  tscale = {}
  plist  = []
  rholist = []
  plots = []

  for i in xrange(len(ts)):
    (cur_plist, cur_rholist) = tscale.get(ts[i],([],[]))
    tscale[ts[i]] = (cur_plist+[ps[i]],cur_rholist+[rhos[i]])

  for t in tscale:
    cur_plist,cur_rholist = tscale[t]
    p, = plt.plot(cur_plist,cur_rholist,linestyle='--',marker='o')
    plots += [(p,t)]
  
  delta_rho = rhos[1]-rhos[0]
  xmin,xmax = plt.xlim()
  ax = plt.gca()
  ax.set_xlim([0.9,1.3])

  plt.plot([xmin,1.3],[0,0], linestyle=':')
  

  lit_ts,lit_rhos,m_ps,lit_ps = np.loadtxt("../output_files/EoS_LJ_6.0.dat",\
                                           unpack = True)

  tscale = {}
  for i in xrange(len(lit_ts)):
    (cur_plist,cur_mlist,cur_rholist) = tscale.get(lit_ts[i],([],[],[]))
    tscale[lit_ts[i]] = (cur_plist+[lit_ps[i]],cur_mlist+[m_ps[i]],\
                         cur_rholist+[lit_rhos[i]])

  for t in tscale:
    cur_plist,cur_mlist,cur_rholist = tscale[t]
    plt.plot(cur_rholist,cur_plist,linestyle='-',color='black')
    plt.plot(cur_rholist,cur_mlist,linestyle='-',color='purple')

  lines, labels = zip(*plots)
  labels = ["T = "+str(labels[x]) for x in xrange(len(labels))]
  plt.legend(lines,labels,loc=2)

  plt.xlabel("$\\rho$")
  plt.ylabel("p")
  
  plt.show()

if __name__=="__main__":main()
