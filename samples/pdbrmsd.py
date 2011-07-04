#! /usr/bin/python

from sys import argv,exit
import numpy as np
from pynumpdb import PDBreader, readPDBFile, readPDBlist, fitCrd
import getopt

def version():
  print """
  pdbrmsd ver. 0.1
  coded by Shiqiao Du in 12,May,2010
  
  ver. 0.2
  revised in 13,May,2010 
  command line interface was renewaled

  ver. 0.3
  revised in 21,Jun,2010
  supported full atom RMSD and other
  """

def usage():
  print """
  usage: <this script> [options] pdb1 pdb2 ... "native.pdb"
      or <this script> [options] -l pdblist "native.pdb" 

  pdblst is a file each of whose line is the path to a pdb file
  e.g.)
  $cat pdblist
  pdb_dir/struct_001.pdb
  pdb_dir/struct_002.pdb
  pdb_dir/struct_003.pdb


  options
  -h, --help: this massage
  -f, --full_atom: calculate full_atom RMSD (no hydrogen)
  -c, --chains: calculate allChain RMSD
  -o, --out=out_file : print output to "out_file"
  """
  exit(2)

version()


# get command line arguments
try:
  opts, args = getopt.getopt(argv[1:],"hcfo:l:", \
              ["help","chains","full_atom","out=","list="])
  # set native structure's file name
  native = args[-1]
except IndexError,getopt.GetoptError:
  usage()

# set other options
allChain = False
onlyCA = True
out_file = "rmsd.log"
pdblist = None

for o,a in opts:
  if o in ("-h","--help"):
    usage()
  if o in ("-c","--chains"):
    allChain = True
  if o in ("-f","--full_atom"):
    onlyCA = False
  if o in ("-o","--out"):
    out_file = a
  if o in ("-l","--list"):
    print "OK"
    pdblist = a

pdbFiles = []
if pdblist:
  fp = file(pdblist)
  for l in fp:
    pdbFiles.append(l.strip())
  fp.close()
else :
  pdbFiles = args[:len(args)-1]

#pdbs,misc,rmsds = readPDBlist(pdbFiles,native,True,allChain,onlyCA)

reader = PDBreader(True,allChain,onlyCA,native)
rmsds = []
for pdb in pdbFiles:
  crd,misc,rmsd = reader.read(pdb,fit=True)
  rmsds.append(rmsd)

result = zip(rmsds,pdbFiles)

#result = []
#native_pdb, misc = readPDBFile(native,)
#for pdb_file in args[:len(args)-1]:
#  pdb,rmsd = fitCrd(native_pdb,readPDBFile(pdb_file)[0])
#  result.append((rmsd, pdb_file))

#result.sort(key = lambda x:float(x[1]))
#result.sort()

fout = file(out_file,"w")
for r in result:
  fout.write("%8.4f    %s\n" % r)
fout.close()
