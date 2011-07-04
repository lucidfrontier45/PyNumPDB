#! /usr/bin/python

from sys import argv, stdout
from os.path import basename
from pynumpdb import getModel,getSeq,threeToOne,getChain

def pdbToFasta(pdb_file,name="STRUCTURE",ID="PDBID", \
                chain="CHAIN",out_file=None):

  pdb_data, ids = getChain(getModel(file(pdb_file).readlines()))
  seq = getSeq(pdb_data)
  seq_one = "".join([threeToOne(s) for s in seq ])

  if out_file:
    out = file(out_file,"w")
  else:
    out = stdout

  if chain == "CHAIN" and ids[0] != " ":
    chain = ids[0]

  out.write(">%s|%s|%s|SEQUENCE\n" % (name,ID,chain))
  for i in range((len(seq_one)/80)+1):
    out.write("%s\n" % seq_one[i*80:(i+1)*80])

if __name__ == "__main__":
  name = "STRUCTURE"
  ID = "PDBID"
  chain = "CHAIN"
  out_file = None
  pdb_file = None
  for i in range(1,len(argv)):
    if argv[i] == "-i" :
      pdb_file = argv[i+1]
      name = basename(pdb_file)[:4].upper()
    if argv[i] == "-n" :
      name = argv[i+1]
    if argv[i] == "-ID" :
      ID = argv[i+1]
    if argv[i] == "-c" :
      chain = argv[i+1]
    if argv[i] == "-o" :
      out_file = argv[i+1]
  pdbToFasta(pdb_file,name,ID,chain,out_file)

