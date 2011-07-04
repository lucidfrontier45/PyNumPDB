#!/usr/bin/python

from sys import argv
from pynumpdb import anm

# set pdbfile name
pdbfile = argv[1]

# create ab instanc of ANM object
model = anm.ANM(pdbfile)

# setup some properties
model.setup()

# perform NMA
model.runNMA()

# save eigen log to "eigen.log"
model.saveEigenLog("eige.log")

# calc B factor
B = model.calcBfactor()

# calc correlation matrix of 1th vibration
cor = self.correlation(1)

# perturbed structure by randomly move along the vibrational directions
traj = model.anmPerturbed(out_name="anm")
