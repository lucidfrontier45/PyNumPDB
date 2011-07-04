#!/usr/bin/python
# Homology modeling by the automodel class
from sys import argv
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
from random import randint

#log.verbose()    # request verbose output
seed = randint(-20000,2) - int(argv[4])
env = environ(seed)  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', "./"+argv[1]]

""" alignment file
>P1;1cm2
structure:1cm2::A::A::::
--MFQQEVTITAPNGLDTRPAAQFVKEAKGFTSEITVTSNGKSASAKSLFKLQTLGLTQGTVVTISAEGEDEQKA
VEHLVKLMAELE----*

>P1;1pfh
structure:1pfh::A::A::::
--MFQQEVTITAPNGL-TRPAAQFVKEAKGFTSEITVTSNGKSASAKSLFKLQTLGLTQGTVVTISAEGEDEQKA
VEHLVKLMAELE----*

>P1;t000
sequence:t000::A::A::::
--MEKKEFHIVAETGIHARPATLLVQTASKFNSDINLEYKGKSVNLKSIMGVMSLGVGQGSDVTITVDGADEAEG
MAAIVETLQKEGLA--*
"""

templates = []
ali_file = argv[2]
with file(ali_file) as f:
  for line in f:
    if line[0] == ">":
      name = line.split(";")[1].strip()
      if name != "t000":
        templates.append(name)
print templates


a = automodel(env,
              alnfile  = ali_file,     # alignment filename
              knowns   = templates,              # codes of the templates
              sequence = 't000')              # code of the target
a.starting_model= int(argv[4])                 # index of the first model
a.ending_model  = a.starting_model - 1 + int(argv[5])                # index of the last model
                                    # (determines how many models to calculate)
a.deviation = float(argv[3])        # max deviation to purtubed initial structure
a.make()                            # do the actual homology modeling

