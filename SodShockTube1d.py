from sph3d import *
from plotter import Plotter1D
from numpy import arange, ones

D = 1

def setSimulationParams():
    Params = ParamsFactory(1)

    S = Params()

    S.time = 0.0
    S.dt = 0.0001
    S.nsteps = 1000
    S.np = 990
    S.eta = 1.2
    S.alpha_min = 0.1
    S.beta = 2.0
    S.epsilon = 0.01
    S.sigma = 0.1
    S.kernel_const = 2./3.
    S.rmax[0] = 1.0
    S.Lboundary[0] = 1
    S.Rboundary[0] = 1

    return S

def setSimulationParticles(S):
    Particle = ParticleFactory(1)

    P = (Particle*S.np)()

    r = arange(3,dtype=float)/2

    dens = ones(2,dtype=float)
    dens[0] *= 8
    dens /= dens.max()

    p = ones(2,dtype=float)
    p[0] *= 10.
    p /= p.max()

    u = eos(dens, p)

    m = densArray1DToParticlePositions(S,P,r,dens)
    uArray1DToParticleU(S,P,r,u)

    for i in range(S.np):
        P[i].m = m
        P[i].h = P[i].m / P[i].rho
        P[i].e = P[i].u+0.5*P[i].v[0]**2
        P[i].Omega = 1.0
        P[i].alpha = S.alpha_min
        P[i].alphau = S.alpha_min
        P[i].valpha = 0.0
        P[i].valphau = 0.0
        P[i].v[0] = 0.0
        P[i].a[0] = 0.0
        P[i].g[0] = 0.0
        P[i].vu = 0.0
        P[i].ve = 0.0
        P[i].adiss[0] = 0.0
        P[i].vudiss = 0.0
        P[i].vediss = 0.0

    return P

def densArray1DToParticlePositions(S,P,r,dens):
    nparticles = (dens/dens.sum() * S.np).astype(int)
    index = 0
    for i in range(dens.size):
        for j in range(nparticles[i]):
            P[index].r[0] = r[i]+(j+1./2)*(r[i+1]-r[i])/nparticles[i]
            P[index].rho = dens[i]
            index += 1
    m = r[1]*dens[0]/nparticles[0]
    
    return m

def uArray1DToParticleU(S,P,r,u):
    for i in range(u.size):
        for j in range(S.np):
            if ((P[j].r[0] >= r[i]) & (P[j].r[0] < r[i+1])):
                P[j].u = u[i]

def eos(rho, P):
    gamma = 7./5.

    return P/((gamma-1)*rho)
