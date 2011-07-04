"""
A package for manipulate sequence
Blast result format

Sbjct: 10  WVVG-KDKPTYDEIFYTLSPVNGKITGANAKKEMVKSKLPNTVLGKIWKLADVDKDGLLD 68
"""

import sys
import re


############################################################################
# sequence manipulation

codes_one = {"G" : "GLY", "A" : "ALA", "L" : "LEU", "I" : "ILE", 
             "R" : "ARG", "K" : "LYS", "M" : "MET", "C" : "CYS",
             "Y" : "TYR", "T" : "THR", "P" : "PRO", "S" : "SER",
             "W" : "TRP", "D" : "ASP", "E" : "GLU", "N" : "ASN",
             "Q" : "GLN", "F" : "PHE", "H" : "HIS", "V" : "VAL",
             "M" : "MSE", "D" : "ASD", "H" : "HID", "H" : "HIP",
             "D" : "GLD"}  

codes_three = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
               "ARG" : "R", "LYS" : "K", "MET" : "M", "CYS" : "C",
               "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
               "TRP" : "W", "ASP" : "D", "GLU" : "E", "ASN" : "N",
               "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V",
               "MSE" : "M", "ASD" : "D", "HID" : "H", "HIP" : "H",
               "GLD" : "D"}

def oneToThree(one_letter):
  o = one_letter.strip().upper()
  if o in codes_one.keys():
    return codes_one[o]
  else:
    sys.stderr.write( "no such code %s\n" % one_letter)
    #return "XXX"
    return ""

def threeToOne(three_letter):
  t = three_letter.strip().upper()
  if t in codes_three.keys():
    return codes_three[t]
  else :
    sys.stderr.write( "no such code %s\n" % three_letter)
    #return "X"
    return ""

def getSeq(pdb_data):
  lastResId = -99999
  seq = []
  for line in pdb_data:
   
    recName = line[0:6].strip()
    if not recName in ("ATOM", "HETATM"):
      continue
  
    resId = int(line[22:26].strip())
    if not resId == lastResId:
      res_name = line[17:20]
      lastResId = resId
      seq.append(res_name)
  
  return seq

