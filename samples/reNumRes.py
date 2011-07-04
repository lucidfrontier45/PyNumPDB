#! /usr/bin/python

from sys import argv
from pynumpdb import reNumRes, getChain, getModel

pdb_data = file(argv[1]).readlines()


re_numbered = reNumRes(getChain(getModel(pdb_data))[0])

for line in re_numbered:
  print line,

