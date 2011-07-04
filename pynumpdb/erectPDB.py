#!/usr/bin/python

import sys
import numpy as np
import pdbtools
import mcmclib

def getVarXY(angles,pdb_crd):
  rotated_crd = pdbtools.rotateCrd(pdb_crd,angles)[:,0:2]
  var = np.cov(rotated_crd.T).trace()
  return var
  
def erectPDB(pdb_file):
  pdb_crd, misc = pdbtools.readPDBFile(pdb_file,True,True)
  init_angle = np.random.rand(3) * 360.0 - 180.0
  print "initial angles = ",
  print init_angle
  ms = mcmclib.MultiMCSystem(getVarXY,init_angle,nsys=10,maxd=10.0, \
                      Tf=20.0,args=(pdb_crd,))
  ms.globalOptimize(max_exchange=100)
  ms.localOptimize()
  return (ms.best_score, ms.best_crd)

if __name__ ==  "__main__":
  result = []
  for i in range(int(sys.argv[2])):
    result.append(erectPDB(sys.argv[1]))
  result.sort(key=lambda x:x[0])
  print """######################################### 
                          result
        """
  for r in result:
    print r
