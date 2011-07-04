#! /usr/bin/python

import re

def getFirstBlastResult(fasta_data):
  score_re = re.compile("^ Score")
  sbj_re = re.compile("^Sbjct:")
  init_count = 0
  sbj_num = []
 
  for line in fasta_data:
    if score_re.search(line):
      init_count += 1
    if init_count == 0:
      continue
    elif init_count  > 1:
      break
    if sbj_re.search(line):
      print line,
      sbj_num.append(line.split()[1])
      sbj_num.append(line.split()[3])

  return sbj_num[0], sbj_num[-1]



if __name__ == "__main__":
  from sys import argv
  num = getFirstBlastResult(file(argv[1]).readlines())
  print num

