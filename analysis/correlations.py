## This code averages the correlations over time. We assume only steady-state measures have been made.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import sys
import glob

dirs = glob.glob('../../output_1DsppBD/%s'%(sys.argv[1]))

for dir_name in dirs:
  filename="%s/corr.dat"%(dir_name)
  if not (os.path.isfile(filename)):
    print("File doesn't exist") 
    continue
  output="%s/avg_corr.dat"%(dir_name)
  #if (os.path.isfile(output)):
  #  print("Average already computed")
  #  continue
  print("Computing average correlations...")
  d=pd.read_csv(filename,sep="\t",header=None,index_col=0)
  C=d.iloc[1,:].copy()
  C=C-C
  
  for i in d.index:
    C=C+d.loc[i,:]/len(d.index)

  C=C[1:len(C)-1]
  C.to_csv(output,sep=" ",index=False,header=False)
