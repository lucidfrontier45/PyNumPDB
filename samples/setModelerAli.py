#! /usr/bin/python


""" clustalw result format
1PTF|PDBID|CHAIN|SEQUENCE      --------------------------------------------------
1MU4|PDBID|CHAIN|SEQUENCE      MVQQKVEVRLKTGLQARPAALFVQEANRFTSDVFLEKDGKKVNAKSIMGL
                                                                                 

1PTF|PDBID|CHAIN|SEQUENCE      ------------------------------------MEKKEFHIVAETGI
1MU4|PDBID|CHAIN|SEQUENCE      MSLAVSTGTEVTLIAQGEDEQEALEKLAAYVQEEVLMVQQKVEVRLKTGL
                                                                   * :::..:  :**:

1PTF|PDBID|CHAIN|SEQUENCE      HARPATLLVQTASKFNSDINLEYKGKSVNLKSIMGVMSLGVGQGSDVTIT
1MU4|PDBID|CHAIN|SEQUENCE      QARPAALFVQEANRFTSDVFLEKDGKKVNAKSIMGLMSLAVSTGTEVTLI
                               :****:*:** *.:*.**: ** .**.** *****:***.*. *::**: 

1PTF|PDBID|CHAIN|SEQUENCE      VDGADEAEGMAAIVETLQKEGLA
1MU4|PDBID|CHAIN|SEQUENCE      AQGEDEQEALEKLAAYVQEEVLQ
                               .:* ** *.:  :.  :*:* * 
"""

"""
C; A sample alignment in the PIR format; used in tutorial

>P1;5fd1
structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA
EVWPNITEKKDPLPDAEDWDGVKGKLQHLER*

>P1;1fdx
sequence:1fdx:1    : :54   : :ferredoxin:Peptococcus aerogenes: 2.00:-1.00
AYVINDSC--IACGACKPECPVNIIQGS--IYAIDADSCIDCGSCASVCPVGAPNPED-----------------
-------------------------------*
"""

from sys import argv, stdout
import re

ali = {}
chain_ids = {}
got_target = False
file_name = argv[1]
try:
  target_id = argv[2]
except:
  target_id = "t000"
out = stdout

fp = file(file_name)
for line in fp:
  if len(line) > 10 and line[4] == "|":
    pdbid = "".join(line[0:4].lower())
    ali.setdefault(pdbid,"")
    ali[pdbid] += (line.strip().split()[1])
    chain_ids.setdefault(pdbid,"")
    curr_chain = line.strip().split("|")[2]
    if curr_chain != "CHAIN":
      chain_ids[pdbid] = curr_chain


fp.close()

for k in ali:
  ali[k] += "*"
  s = re.sub("\.","-",ali[k].upper())
  c = chain_ids[k]
  out.write(">P1;%s\n" % k)
  if k == target_id:
    seq_type = "sequence"
  else :
    seq_type = "structure"
  out.write("%s:%s::%s::%s::::\n" % (seq_type, k, c, c ))
  for i in range((len(s)/75)+1):
    out.write("%s\n" % s[i*75:(i+1)*75])
  out.write("\n")

