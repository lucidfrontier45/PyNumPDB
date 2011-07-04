#!/usr/bin/python

from sys import argv,exit,stdout
import getopt
import numpy as np
import pynumpdb
import pynumpdb.clustering as pdbcl

def version():
  pdbcl.version()

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
  -c, --clusters=n : set number of clusters, n must be an integer
  -t, --tops=n : retrieve the n structures nearest to largest cluster centroid
                 n must be an integer
  -r, --ref=ref_pdb : set reference PDB
  -m, --mode=method : kemans and gmm are supported
  -n. --norm=method : specify normalization method i.e. "pca", "whiten" and "pca_whiten"
  """
  exit(2)

def main():
  

  version()

  # get command line arguments
  try:
    opts, args = getopt.getopt(argv[1:],"hc:t:o:r:m:n:i:d:p:", \
            ["help","clusters=","tops=","out=","ref=","mode=", \
             "norm=","iter=","data=","pca="])
  except getopt.GetoptError:
    usage()

  
  # set other options
  nc = 4
  tops = 5
  out_file = "clustlog"
  refPDB_file = None
  mode = "kmeans"
  norm = "pca"
  norm_dim = 20
  itermax = 500
  data_file = None

  for o,a in opts:
    if o in ("-h","--help"):
      usage()
    if o in ("-c","--clusters"):
      nc = int(a)
    if o in ("-t","--tops"):
      tops = int(a)
    if o in ("-o","--out"):
      out_file = a
    if o in ("-r","--ref"):
      refPDB_file = a
    if o in ("-n","--norm"):
      norm = a
    if o in ("-m","--mode"):
      mode = a
    if o in ("-i","--iter"):
      itermax = int(a)
    if o in ("-d","--data"):
      data_file = a
    if o in ("-p","--pca"):
      norm_dim = int(a)

  if len(args) < 1:
    usage()
  else:
    pdb_list = args[0]
    pdbfiles = file(pdb_list).readlines()
    
  # set input file
  if data_file:
    # read  from data_file
    data = np.loadtxt(data_file)
  else:
    # read pdblist file, fit each structure to refPDB and normalize
    pdbdata, miscs, rmsds = pynumpdb.readPDBlist(pdbfiles,refPDB_file)
    data = pdbcl.normalize(pdbdata,norm,norm_dim)
    np.savetxt(out_file+"_data",data)

  # run main routines
  cm = pdbcl.ClusteringMachine(nc)
  try:
    res = cm.clustering(data,mode,itermax)
  except ClusteringModeError, m:
    print m
    exit(2)

  if mode == "gmm":
    #log_likelihood = res[2][-1]
    log_likelihood = res[1].sum()
    comp = pdbcl.complexity(data,nc)
    aic = (-log_likelihood + comp) * 2.0
    bic = -log_likelihood * 2.0 + comp * np.log(len(pdbfiles))
    print "nc = %3d, aic = %8.3e, bic = %8.3e" % (nc,aic,bic)

  new_codes, pop = cm.rankCodes(res[0])
  clusts = cm.sortedCluster(new_codes,res[1])

  # print sequences
  out = file(out_file+".seq","w")
  out.write("T= %d\n" %len(new_codes))
  for i in new_codes:
    out.write("%d " % (i+1))
  out.write("\n")
  out.close()

  # print detail
  for i in range(len(pop)):
    out = file("%s_%03d" %(out_file,i),"w")
    out.write("# population of clust %d is %8.3f\n" %(i,pop[i]))
    for c in clusts[i]:
      out.write("%3d %8.3f %8.3f %8.3f %8.3f %s\n" %(i,c[1], \
        data[c[0],0], data[c[0],1], data[c[0],2], pdbfiles[c[0]].strip()))
    out.write("\n")
    out.close()

  # generate plot script
  out = file("%s.plt" %(out_file),"w")
  out.write(
"""#!/usr/bin/gnuplot
set xlabel "PC1"
set ylabel "PC2"
#set term postscript color enhanced 20
plot""" 
)
  for i in range(len(clusts)):
    out.write("'%s_%03d' u 3:4 title 'clust %03d'" %(out_file,i,i))
    if i < len(clusts) - 1:
      out.write(", ")
  out.close()

if __name__ == "__main__":
  main()
