## This code computes J and M order parameters from sizedistr.dat and nclusters.dat (average).

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import sys
import glob

dirs = glob.glob('../../output_1DsppBD/%s'%(sys.argv[1]))

OP=pd.DataFrame(data=np.zeros((len(dirs),5)));
i=0

##../../output_1DsppBD/sim_a0.001_f0.800_t0010000000_L00500_D1.000_Fp0.00_eps0.04167_CMOB1_IS0_tint10000

for dir_name in dirs:
  alpha=float(dir_name[27:31])
  phi=float(dir_name[34:38])
  T=int(dir_name[41:50])
  L=int(dir_name[53:57])
  N=L*phi
  D=float(dir_name[60:64])
  beta=float(pd.read_csv("%s/parameters.dat"%(dir_name),engine='python',sep=": ",header=None,names=["Par","Value"])["Value"][11])
  v=float(pd.read_csv("%s/parameters.dat"%(dir_name),engine='python',sep=": ",header=None,names=["Par","Value"])["Value"][10])
  ##v=float(dir_name[68:71])
  eps=float(dir_name[76:82])
  CMOB=int(dir_name[87])
  IS=int(dir_name[91])
  Tint=int(dir_name[97:102])
  file_sizes="%s/sizedistr.dat"%(dir_name)
  file_nclust="%s/nclusters.dat"%(dir_name)
  if not (os.path.isfile(file_sizes)):
    print("File doesn't exist") 
    continue
  if not (os.path.isfile(file_nclust)):
    print("File doesn't exist") 
    continue
  if i==0:
    output="../../output_1DsppBD/OP_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(T,L,v,eps,CMOB,IS,Tint)
  print("Computing order parameters...")
  d_size=pd.read_csv(file_sizes,sep="\t",header=None,names=["l","p"])
  d_clust=pd.read_csv(file_nclust,sep="\t",header=None,index_col=0)
  Nc=d_clust.loc[Tint:,1]
  Nc_avg=sum(Nc)/len(Nc)
  p=d_size["p"]
  l=d_size["l"]
  n=p*Nc_avg
  pl=p*l
  nl=n*l
  J_this=sum(nl)/N
  M_this=sum(pl)/N

  OP.iloc[i,0]=alpha
  OP.iloc[i,1]=phi
  OP.iloc[i,2]=beta
  OP.iloc[i,3]=J_this
  OP.iloc[i,4]=M_this
  
  i=i+1


OP=OP[OP[0]>0]
JM=OP.copy()

alpha=[]
for a in OP[0]:
  if not a in alpha:
    alpha.append(a)
alp=np.array(alpha.copy())
alpha=np.array(alpha)
alpha=alpha-alpha
for i in range(len(alpha)):
  alpha[i]=min(alp[alp>max(alpha)])

phi=[]
for f in OP[1]:
  if not f in phi:
    phi.append(f)
ph=np.array(phi.copy())
phi=np.array(phi)
phi=phi-phi
for i in range(len(phi)):
  phi[i]=min(ph[ph>max(phi)])

beta=[]
for b in OP[2]:
  if not b in beta:
    beta.append(b)
be=np.array(beta.copy())
beta=np.array(beta)
beta=beta-beta
for i in range(len(beta)):
  beta[i]=min(be[be>max(beta)])

print("Fp=%.2f\neps=%.5f\n"%(v,eps))

print("\n\nConstant beta\n")

#CONSTANT BETA:
for b in beta:
  output_J=open("../../output_1DsppBD/J_constBeta_D%.5f_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(1/b,T,L,v,eps,CMOB,IS,Tint), "w")
  output_M=open("../../output_1DsppBD/M_constBeta_D%.5f_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(1/b,T,L,v,eps,CMOB,IS,Tint), "w")
  for a in alpha:
    print("Tumbling rate: %f"%(a))
    if any(OP[OP[2]==b][0]==a):
      for f in phi:
        print("  Density: %f"%(f))
        if any(OP[OP[2]==b][OP[0]==a][1]==f):
          print("    We have value")
          output_J.write("%f %f %f\n"%(a,f,OP[OP[2]==b][OP[0]==a][OP[1]==f][3]))
          output_M.write("%f %f %f\n"%(a,f,OP[OP[2]==b][OP[0]==a][OP[1]==f][4]))
        else:
          print("    We don't have value")
          output_J.write("%f %f NaN\n"%(a,f))
          output_M.write("%f %f NaN\n"%(a,f))
    else:
      for f in phi:
        print("  Density: %f"%(f))
        print("    We don't have value")
        output_J.write("%f %f NaN\n"%(a,f))
        output_M.write("%f %f NaN\n"%(a,f))
    output_J.write("\n")
    output_M.write("\n")
  output_J.close()
  output_M.close()

print("\n\nConstant alpha\n")

#CONSTANT ALPHA:
for a in alpha:
  output_J=open("../../output_1DsppBD/J_constAlpha_a%.5f_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(a,T,L,v,eps,CMOB,IS,Tint), "w")
  output_M=open("../../output_1DsppBD/M_constAlpha_a%.5f_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(a,T,L,v,eps,CMOB,IS,Tint), "w")
  for f in phi:
    if any(OP[OP[0]==a][1]==f):
      for b in beta:
        if any(OP[OP[0]==a][OP[1]==f][2]==b):
          output_J.write("%f %f %f\n"%(f,b,OP[OP[0]==a][OP[1]==f][OP[2]==b][3]))
          output_M.write("%f %f %f\n"%(f,b,OP[OP[0]==a][OP[1]==f][OP[2]==b][4]))
        else:
          output_J.write("%f %f NaN\n"%(f,b))
          output_M.write("%f %f NaN\n"%(f,b))
    else:
      for f in phi:
        output_J.write("%f %f NaN\n"%(f,b))
        output_M.write("%f %f NaN\n"%(f,b))
    output_J.write("\n")
    output_M.write("\n")
  output_J.close()
  output_M.close()

print("\n\nConstant phi\n")

#CONSTANT PHI:
for f in phi:
  output_J=open("../../output_1DsppBD/J_constPhi_f%.5f_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(f,T,L,v,eps,CMOB,IS,Tint), "w")
  output_M=open("../../output_1DsppBD/M_constPhi_f%.5f_t%.10d_L%.5d_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d.dat"%(f,T,L,v,eps,CMOB,IS,Tint), "w")
  for a in alpha:
    if any(OP[OP[1]==f][0]==a):
      for b in beta:
        if any(OP[OP[1]==f][OP[0]==a][2]==b):
          output_J.write("%f %f %f\n"%(a,b,OP[OP[1]==f][OP[0]==a][OP[2]==b][3]))
          output_M.write("%f %f %f\n"%(a,b,OP[OP[1]==f][OP[0]==a][OP[2]==b][4]))
        else:
          output_J.write("%f %f NaN\n"%(a,b))
          output_M.write("%f %f NaN\n"%(a,b))
    else:
      for f in phi:
        output_J.write("%f %f NaN\n"%(a,b))
        output_M.write("%f %f NaN\n"%(a,b))
    output_J.write("\n")
    output_M.write("\n")
  output_J.close()
  output_M.close()
