#! /usr/bin/python

from sys import argv
from pynumpdb import getPeptide


pdb_file = argv[1]

if argv[2] == "-ali" :
  blast_data = file(argv[3]).readlines()
  start, end = getFirstBlastResult(blast_data) 
  start = int(start)
  end = int(end)
else:
  start = int(argv[2])
  end = int(argv[3])

reReNum = False
noh = True

if "-r" in argv:
  reReNum = True
if "-h" in argv:
  noh = False
  
peptides = getPeptide(file(pdb_file).readlines(),start,end,reReNum,noh)
for p in peptides:
  print p,
