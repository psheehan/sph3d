from sph3d import *
from numpy import arange, meshgrid, ones, pi
from numpy.random import random

def setSimulationParams():
    Params = ParamsFactory(2)

    S = Params()

    S.time = 0.0
    S.dt = 0.01
    S.nsteps = 100
    S.np = 1089
    S.eta = 1.2
    S.alpha_min = 0.1
    S.beta = 2.0
    S.epsilon = 0.01
    S.sigma = 0.1
    S.kernel_const = 10./(7*pi)
    S.rmax[0] = 1.0
    S.rmax[1] = 1.0
    S.Lboundary[0] = 1
    S.Lboundary[1] = 1
    S.Rboundary[0] = 1
    S.Rboundary[1] = 1

    return S

def setSimulationParticles(S):
    Particle = ParticleFactory(2)

    P = (Particle*S.np)()

    nx = 1
    ny = 1

    x = arange(nx+1,dtype=float)/nx
    y = arange(ny+1,dtype=float)/ny

    dens = ones((nx,ny),dtype=float)

    u = ones((nx,ny),dtype=float)*1.0e-5

    m = densArray2DToParticlePositions(S,P,x,y,dens,nx,ny)
    uArray2DToParticleU(S,P,x,y,u,nx,ny)

    for i in range(S.np):
        if (((P[i].r[0]-0.5)**2+(P[i].r[1]-0.5)**2)**(1./2) <= 0.05):
            P[i].u = 1.0
        P[i].m = m
        P[i].h = (P[i].m / P[i].rho)**(1./2)
        P[i].Omega = 1.0
        P[i].alpha = S.alpha_min
        P[i].alphau = S.alpha_min
        P[i].valpha = 0.0
        P[i].valphau = 0.0
        P[i].v[0] = 0.0
        P[i].v[1] = 0.0
        P[i].a[0] = 0.0
        P[i].a[1] = 0.0
        P[i].g[0] = 0.0
        P[i].g[1] = 0.0
        P[i].e = P[i].u+0.5*(P[i].v[0]**2+P[i].v[1]**2)
        P[i].vu = 0.0
        P[i].ve = 0.0
        P[i].adiss[0] = 0.0
        P[i].adiss[1] = 0.0
        P[i].vudiss = 0.0
        P[i].vediss = 0.0

    return P

def densArray2DToParticlePositions(S,P,x,y,dens,nx,ny):
    nparticles = (dens/dens.sum() * S.np).astype(int)
    nparticles_x = (nparticles**(1./2)).astype(int)
    nparticles_y = (nparticles**(1./2)).astype(int)

    index = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nparticles_x[i,j]):
                for l in range(nparticles_y[i,j]):
                    if (l%2 == 0):
                        P[index].r[0] = x[i]+(k+1./2)*(x[i+1]-x[i])/ \
                            nparticles_x[i,j]
                    else:
                        P[index].r[0] = x[i]+k*(x[i+1]-x[i])/ \
                            nparticles_x[i,j]
                    P[index].r[1] = y[j]+(l+1./2)*(y[j+1]-y[j])/ \
                        nparticles_y[i,j]
                    P[index].rho = dens[i,j]
                    index += 1
    m = x[1]*y[1]*dens[0,0]/nparticles[0,0]
    
    return m

def uArray2DToParticleU(S,P,x,y,u,nx,ny):
    for i in range(nx):
        for j in range(ny):
            for k in range(S.np):
                if (((P[k].r[0] >= x[i]) & (P[k].r[0] < x[i+1])) &
                    ((P[k].r[1] >= y[j]) & (P[k].r[1] < y[j+1]))):
                    P[k].u = u[i,j]

def eos(rho, P):
    gamma = 7./5.

    return P/((gamma-1)*rho)
