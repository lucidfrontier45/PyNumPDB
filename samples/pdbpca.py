#!/usr/bin/python

from sys import argv,exit
import getopt
import numpy as np
import pynumpdb
import pynumpdb._pca
from scikits.learn.decomposition import PCA

def version():
  print """
  pdbpca ver. 0.2
  coded by Shiqiao Du in 24,April,2010
  revised in 4,May,2010
  """

def usage():
  print """
  usage: <this script> [options] pdblist
  
  pdblst is a file each of whose line is the path to a pdb file
  e.g.)
  $cat pdblist
  pdb_dir/struct_001.pdb
  pdb_dir/struct_002.pdb
  pdb_dir/struct_003.pdb


  options
  -h, --help: this massage
  -o, --out=out_file : print output to "out_file"
  -r, --ref=ref_pdb : set reference PDB
  -n. --npcs=n : compute and print n principal components, n must be an integer
  """
  exit(2)

def PDBpca(pdblist_file, npcs=5,refPDB_file=None):

  # read pdblist file and fit each structure to refPDB
  pdbdata, miscs, rmsds = pynumpdb.readPDBlist(pdblist_file,refPDB_file)
   
  # run PCA
  #v,P,PC = pynumpdb._pca.pca_train(pdbdata,npcs,do_norm=0)
  pca = PCA()
  pca.fit(pdbdata)
  v = pca.explained_variance_
  P = pca.components_
  PC = pca.transform(pdbdata)
  print v
  print P
  print len(PC),len(PC[0])
  #print PC.T

  return v,P,PC


def main():
  

  version()

  # get command line arguments
  try:
    opts, args = getopt.getopt(argv[1:],"ho:r:n:", \
                ["help","out=","ref=","npcs="])
  except getopt.GetoptError:
    usage()

  # set input pdb list
  if len(args) < 1:
    usage()
  else:
    pdb_list = file(args[0]).readlines()

  
  # set other options
  out_file = "pca.log"
  npcs = 5
  refPDB_file = None

  for o,a in opts:
    if o in ("-o","--out"):
      out_file = a
    if o in ("-n","--npcs"):
      npcs = int(a)
    if o in ("-r","--ref"):
      refPDB_file = a

  # run main routine
  v,P,PC = PDBpca(pdb_list,npcs,refPDB_file)

  # save PCs to out_file
  print "saveing result to %s" % out_file
  np.savetxt(out_file,PC,fmt="%6.3f")

if __name__ == "__main__":
  main()
