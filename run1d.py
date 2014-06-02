#!/usr/bin/env python3

from sph3d import *
from plotter import Plotter1D
from time import time

# Import the simulation to be run.

from SodShockTube1d import *

# Now do all of the important stuff.

S = setSimulationParams()

P = setSimulationParticles(S)

SPH = SPHFactory(1)
H = SPH(S,P)

H.initialize()

Plot = Plotter1D(S,P)

Plot.plot()

for i in range(S.nsteps):
    H.step()
    Plot.plot()

Var = raw_input("type enter to quit")
