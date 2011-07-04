import numpy as np
from scipy import linalg
from _selection import getModel, getChain
from sys import stdout

"""
A package for IO and numerical operation on PDB
"""

class PDBreader:
  """
  class for reading PDB files
  """
  def __init__(self,center=True,allChain=False,onlyCA=True,refPDB=None):
    """
    args
    center : boolean, if perform centering
    allChain : boolean, if use all chains
    onlyCA : boolean, if use only CA
    refPDB : strings, path to reference PDB file
    """
    self.center = center
    self.allChain = allChain
    self.onlyCA  = onlyCA
    self.refPDB = refPDB
    self.ref_crd = None
    self.ref_misc = None
    if self.refPDB:
      self.ref_crd, self.ref_misc = readPDBFile( \
          self.refPDB,self.center,self.allChain,self.onlyCA)
    
  def read(self,pdbfile,fit=False):
    rmsd = 0.0
    crd, misc = readPDBFile(pdbfile,self.center,self.allChain,self.onlyCA)
    if fit:
      crd, rmsd = fitCrd(self.ref_crd,crd)
    return crd,misc,rmsd

  def readList(self,pdblist,fit=True):
    pdbfiles = file(pdblist).readlines()
    crds,miscs,rmsds = readPDBlist(pdbfiles,self.ref_crd,fit, \
        self.allChain,self.onlyCA)
    return crds, miscs, rmsds


def centeringCrd(crd):
  temp = np.array(crd)
  return temp - temp.mean(0)

def setRotateMatrix(angle, axis="z"):
  r = np.zeros(9).reshape(3,3)
  rad = angle / np.pi
  if axis in ("x","X"):
    r[0,0] = np.cos(angle)
    r[0,1] = np.sin(angle)
    r[1,0] = -np.sin(angle)
    r[1,1] = np.cos(angle)
    r[2,2] = 1.0
  elif axis in ("y","Y"):
    r[0,0] = np.cos(angle)
    r[0,2] = np.sin(angle)
    r[2,0] = -np.sin(angle)
    r[2,2] = np.cos(angle)
    r[1,1] = 1.0
  elif axis in ("z","Z"):  
    r[1,1] = np.cos(angle)
    r[1,2] = np.sin(angle)
    r[2,1] = -np.sin(angle)
    r[2,2] = np.cos(angle)
    r[0,0] = 1.0
  else:
    r = np.identity(3)

  return r

def rotateCrd(crd,angles):
  B = setRotateMatrix(angles[0],"x")
  C = setRotateMatrix(angles[1],"z")
  D = setRotateMatrix(angles[2],"x")
  A = np.dot(B,np.dot(C,D))
  return np.dot(crd,A.T)

def pdbToArray(data,center=True,allChain=False,onlyCA=True,noh=True):
  crd = []
  misc = []
  pdb_data = getModel(data,noh=noh)
  if allChain:
    print("use allChain\n") 
  else:
    pdb_data,ids = getChain(pdb_data)
  for line in pdb_data:
    #if line[0:3] == "TER" : break
    if line[12:16].strip() == "CA" or not onlyCA:
      try:
        crd.append([np.double(line[30:38]),
                    np.double(line[38:46]),
                    np.double(line[46:54])])
        misc.append((line[:30],line[54:]))
      except:
        exit()
  if center:
    crd = centeringCrd(crd)
  print "length = %d" % len(crd)
  return np.array(crd),misc

def convCrd(crd,n=1):
  """ convert N 3-dim crds to 1 3N-dim crd and vice versa""" 
  if(n == 1):
    return crd.flatten()
  elif(n == 3):
    return crd.reshape(crd.size/3,3)
  else:
    print "warning invalid dimension %d" %n
    return crd

def fitCrd(crd1,crd2):
  ccrd1 = convCrd(crd1,3)
  ccrd2 = convCrd(crd2,3)

  #fit crd2 to crd1 by SVD
  R = np.dot(ccrd2.T,ccrd1)
  Wt, S, V = linalg.svd(R)
  U = np.dot(V.T,Wt.T)

  # calc RMSD
  rmsd = np.sqrt(np.abs(\
            np.dot(ccrd1.flatten(),ccrd1.flatten()) + \
            np.dot(ccrd2.flatten(),ccrd2.flatten()) - \
            2.0 * S.sum())/len(ccrd1))
  return np.dot(ccrd2,U.T),rmsd

def readPDBFile(file_name,center=True,allChain=False,onlyCA=True,noh=True):
  print "readind " + file_name, 
  pdb_data = file(file_name).readlines()
  crd,misc = pdbToArray(pdb_data,center,allChain,onlyCA,noh)
  return convCrd(crd,1),misc

def readPDBlist(pdblist,refPDB_file=None,if_fit=True, \
                allChain=False,onlyCA=True,noh=True):
  pdbtemp = []
  miscs = []
  rmsds = []
 
  # get reference PDB
  if refPDB_file:
    refcrd, misc = readPDBFile(refPDB_file.strip(),if_fit,allChain,onlyCA,noh)
  else :
    pdbfile = pdblist[0]
    refcrd, misc = readPDBFile(pdbfile.strip(),if_fit,allChain,onlyCA,noh)
  reflen = len(refcrd)

  # loop over pdblists 
  for pdbfile in pdblist:
    crd, misc = readPDBFile(pdbfile.strip(),if_fit,allChain,onlyCA,noh)
    if not len(crd) == reflen:
      print "length not match, skip this file"
      continue
    # fit to refPDB and append to temp array
    miscs.append(misc)
    if if_fit:
      fitted_crd,rmsd = fitCrd(refcrd,crd)
      rmsds.append(rmsd)
      pdbtemp.append(convCrd(fitted_crd,1))
    else :
      pdbtemp.append(convCrd(crd,1))

  pdbdata = np.array(pdbtemp)
  
  return pdbdata,misc,rmsds

def writePDBFile(file_name,crd,misc):
  ccrd = convCrd(crd,3)
  fp = file(file_name,"w")
  for i in range(len(misc)):
    fp.write("%s%8.3f%8.3f%8.3f%s" % (misc[i][0], ccrd[i][0], \
      ccrd[i][1], ccrd[i][2], misc[i][1]))
  fp.write("TER\n")
  fp.close()

class CAModel():
  def __init__(self,pdb_file,allChain=True):
    self.pdb_file = pdb_file
    self.crd, self.misc = readPDBFile(pdb_file,allChain)
    self.crd = convCrd(self.crd,3)

  def __len__(self):
    return len(self.crd)
