import matplotlib.pyplot as plt
from numpy import arange, concatenate, linspace, log10, zeros
from scipy.interpolate import griddata

class Plotter1D:

    def __init__(self, S, P):
        self.S = S
        self.P = P

        self.f, self.ax = plt.subplots(4, figsize=(7,10), \
                sharex=True, sharey=False)
        self.line = [1,2,3,4]
        plt.ion()

        self.line[0], = self.ax[0].plot([], [], 'bo', markeredgecolor='b')
        self.ax[0].set_xlim(0.0, 1.0)
        self.ax[0].set_ylabel(r'$P$', fontsize = 16)
        self.ax[0].set_ylim(0.0, 1.1)

        self.line[1], = self.ax[1].plot([], [], 'ro', markeredgecolor='r')
        self.ax[1].set_ylabel(r'$\rho$', fontsize = 16)
        self.ax[1].set_ylim(0.0, 1.1)

        self.line[2], = self.ax[2].plot([], [], 'go', markeredgecolor='g')
        self.ax[2].set_ylabel(r'$u$', fontsize = 16)
        self.ax[2].set_ylim(1.5, 3.0)

        self.line[3], = self.ax[3].plot([], [], 'yo', markeredgecolor='y')
        self.ax[3].set_xlabel(r'$r$', fontsize = 16)
        self.ax[3].set_ylabel(r'$v$', fontsize = 16)
        self.ax[3].set_ylim(0.0,1.1)

        self.f.subplots_adjust(hspace=0.2)

        bbprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        self.text = self.ax[0].text(0.75, 1.25, "", \
                transform=self.ax[0].transAxes, fontsize=12, style='normal', \
                bbox=bbprops)

    def dataAsArrays(self):
        r = zeros(self.S.np)
        P = zeros(self.S.np)
        rho = zeros(self.S.np)
        u = zeros(self.S.np)
        v = zeros(self.S.np)

        for i in range(self.S.np):
            r[i] = self.P[i].r[0]
            P[i] = self.P[i].P
            rho[i] = self.P[i].rho
            u[i] = self.P[i].u
            v[i] = self.P[i].v[0]

        return r, P, rho, u, v

    def plot(self):
        r, P, rho, u, v = self.dataAsArrays()

        self.line[0].set_data(r, P)
        self.line[1].set_data(r, rho)
        self.line[2].set_data(r, u)
        self.line[3].set_data(r, v)

        textstr = "$t = %6.4f$\n$dt = %7.2e$"%(self.S.time,self.S.dt)
        self.text.set_text(textstr)

        plt.show()
        plt.draw()

class Plotter2D:

    def __init__(self,S,P):
        self.S = S
        self.P = P

        self.f, self.ax = plt.subplots(nrows=1, ncols=1, figsize=(9,7), \
                sharex=False, sharey=False)
        plt.ion()

        self.ax.set_xlabel(r'$x$', fontsize = 16)
        self.ax.set_ylabel(r'$y$', fontsize = 16)
        self.ax.set_xlim(0.0,1.0)
        self.ax.set_ylim(0.0,1.0)

        bbprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        self.text = self.ax.text(0.75, 1.05, "", \
                transform=self.ax.transAxes, fontsize=12, style='normal', \
                bbox=bbprops)

    def dataAsArrays(self):
        x = zeros(self.S.np)
        y = zeros(self.S.np)
        r = zeros(self.S.np)
        rho = zeros(self.S.np)
        P = zeros(self.S.np)

        for i in range(self.S.np):
            x[i] = self.P[i].r[0]
            y[i] = self.P[i].r[1]
            rho[i] = self.P[i].rho
            P[i] = self.P[i].P

        r = ((x-0.5)**2+(y-0.5)**2)**(1./2)

        return x, y, log10(P)

    def init(self):
        x, y, P = self.dataAsArrays()

        self.scat = self.ax.scatter(x, y, c=P, s=50, edgecolors='none', \
                vmin = -6, vmax=1)
        plt.colorbar(self.scat)

        textstr = "$t = %6.4f$\n$dt = %7.2e$"%(self.S.time,self.S.dt)
        self.text.set_text(textstr)

        plt.show()
        plt.draw()

    def plot(self):
        x, y, P= self.dataAsArrays()

        self.scat.set_offsets(concatenate((x.reshape(x.size,1), \
                y.reshape(y.size,1)),axis=1))
        self.scat.set_array(P)

        textstr = "$t = %6.4f$\n$dt = %7.2e$"%(self.S.time,self.S.dt)
        self.text.set_text(textstr)

        plt.draw()
