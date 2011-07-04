import numpy as np
from scipy import linalg

def canonicalNormalize(data):
  mean = np.mean(data,axis=0)
  stdev = np.std(data,axis=0)
  norm = (data - mean ) / stdev
  return norm

def pca_train(input_data,npcs=5,do_norm=1,do_project=1):
  if do_norm :
    print "normalising data",
    data = canonicalNormalize(input_data)
    print "done" 
  else :
    data = input_data

  print "calculating cov matrix",
  COV = -np.cov(data.T)
  print "done" 
  print len(COV),len(COV[0])

  print "diagnalizing",
  v,P = linalg.eigh(COV,eigvals=(1,npcs))
  print "done" 
  if do_project :
    return -v,P,np.dot(data,P)
  else :
    return -v,P

def pca_proj(P,input_data,do_norm=1,):
  if do_norm :
    data = canonicalNormalize(input_data)
  else :
    data = input_data
  return np.dot(data,Po)

