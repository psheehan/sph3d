#!/usr/bin/env python

from sph3d import *
from plotter import Plotter2D

# Import the simulation to be run.

#from SodShockTube2d import *
from SedovTaylor2d import *
#from KelvinHelmholtz2d import *

# Now do all of the important stuff.

S = setSimulationParams()

P = setSimulationParticles(S)

SPH = SPHFactory(2)
H = SPH(S,P)

H.initialize()

Plot = Plotter2D(S,P)

Plot.init()

for i in range(S.nsteps):
    H.step()
    Plot.plot()

Var = raw_input("type enter to quit")
