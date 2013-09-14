#ifndef GRID_CC
#define GRID_CC

#include <cmath>
#include <vector>
#include "particle.cc"

struct Cell {
    int nparticles;
    int nsources;
    std::vector <int> particle_keys;
};

struct Grid {
    int ncells;
    Cell *c;

    template <int D>
    Grid(int _ncells, Particle<D> *P, int np) {
        ncells = 1;
        for (int i=0; i<D; i++)
            ncells *= _ncells;

        c = new Cell[ncells]();
        for (int i=0; i<ncells; i++)
            c[i].nparticles = 0;

        for (int i=0; i<np; i++) {
            int index[D];
            int key = 0;
            for (int j=0; j<D; j++) {
                index[j] = int(P[i].r[j] * _ncells);
                key += index[j]*pow(_ncells,D-1-j);
            }
            c[key].particle_keys.push_back(i);
            c[key].nparticles++;
        }
    }

    template <int D>
    Particle<D> *get_sinks(int index, Particle<D> *particles) {
        Particle<D> *sinks = new Particle<D>[c[index].nparticles];

        for (int i=0; i<c[index].nparticles; i++)
            sinks[i] = particles[c[index].particle_keys[i]];

        return sinks;
    }

    template <int D>
    void update_sinks(int index, Particle<D> *particles, Particle<D> *sinks) {
        for (int i=0; i<c[index].nparticles; i++)
            particles[c[index].particle_keys[i]] = sinks[i];

        delete sinks;
    }

    template <int D>
    Particle<D> *get_sources(int index, Particle<D> *particles) {
        double h_max = 0.0;
        for (int i=0; i<c[index].nparticles; i++) {
            if (particles[c[index].particle_keys[i]].h > h_max)
                h_max = particles[c[index].particle_keys[i]].h;
        }

        int nacross = 1+int(2*h_max*ncells);

        int start = index - nacross;
        int fin = index + nacross;
        int nsources = 0;

        for (int i=start; i<=fin; i++) {
            if (i >= ncells)
                nsources += c[i-ncells].nparticles;
            else if (i < 0)
                nsources += c[i+ncells].nparticles;
            else
                nsources += c[i].nparticles;
        }
        c[index].nsources = nsources;

        Particle<D> *sources = new Particle<D>[nsources];

        int ind = 0;
        for (int i=start; i<=fin; i++) {
            int cell_index = i;
            if (i >= ncells)
                cell_index -= ncells;
            else if (i < 0)
                cell_index += ncells;

            for (int j=0; j<c[cell_index].nparticles; j++) {
                sources[ind] = particles[c[cell_index].particle_keys[j]];
                ind++;
            }
        }

        return sources;
    }
};

#endif
