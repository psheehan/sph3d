#ifndef PARTICLE_CC
#define PARTICLE_CC

#include "vector.cc"

template<int D>
struct Particle {
    double rho, m;
    Vector<double, D> r;
    Vector<double, D> v;
    Vector<double, D> a;
    Vector<double, D> g;
    double vu;
    double ve;
    Vector<double, D> adiss;
    double vudiss;
    double vediss;
    double h;
    double P;
    double u;
    double e;
    double cs;
    double alpha;
    double alphau;
    double valpha;
    double valphau;
    double Omega;
};

#endif
