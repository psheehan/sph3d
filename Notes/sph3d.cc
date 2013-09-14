#include <cmath>

#define VS(a,s) for (int d=0; d<D; d++) a[d] s

#define VV(a,b) for (int d=0; d<D; d++) a[d] b[d]

#define VSV(a,s,b) for (int d=0; d<D; d++) a[d] s b[d]

#define VVV(a,b,c) for (int d=0; d<D; d++) a[d] b[d] c[d]

#define D 3
#if D==1
# define NORM(a) fabs(a)
#elif D==2
# define NORM(a) sqrt(a[0]*a[0]+a[1]*a[1])
#elif D==3
# define NORM(a) sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#endif


typedef struct {
    double rho, m;
    double r[D];
    double v[D];
    double a[D];
    double g[D];
    double h;
    double P;
    double u;
    double Omega;
} Particle;

typedef struct {
    double dt;
    double eta;
    int np;
} Params;

class SPH {
public:
    Particle *P;
    Params *S;

    SPH(Params *_S, Particle *_P): S(_S), P(_P) {};

    double M_4(double q);
    double dM_4(double q);
    double kernel(double r, double h);
    double grad_kernel(double r, double h);
    double eos(double u, double rho);
    void calculate_P();
    void calculate_rho(Particle *sinks, int nsinks, Particle *sources, 
            int nsources);
    void calculate_accel(Particle *sinks, int nsinks, Particle *sources,
            int nsources);
    void drift(double dt);
    void kick(double dt);
    void update_timestep();
    void step();
};

void SPH::calculate_P() {
    for (int i=0; i<S->np; i++) {
        P[i].P = eos(P[i].u, P[i].rho);
    }
};

double SPH::M_4(double q) {
    if ((0 <= q) && (q < 1)) {
        return 1./4.*(2.-q)*(2.-q)*(2.-q)-(1.-q)*(1.-q)*(1.-q);
    }
    else if ((1 <= q) && (q <= 2)) {
        return 1./4.*(2-q)*(2-q)*(2-q);
    }
    else {
        return 0.0;
    }
};

double SPH::kernel(double r, double h) {

    double sigma;

    if (D == 1) {
        sigma = 2./3.;
    }
    else if (D == 2) {
        sigma = 10./(7.*3.1415927);
    }
    else if (D == 3) {
        sigma = 1./3.1415927;
    }

    return 1./pow(h,D)*M_4(r/h)*sigma;
};

double SPH::dM_4(double q) {
    if ((0 <= q) && (q < 1)) {
        return -3.+9./4.*q;
    }
    else if ((1 <= q) && (q <= 2)) {
        return -3./4.*(2.-q)*(2.-q)/q;
    }
    else {
        return 0.0;
    }
};

double SPH::grad_kernel(double r, double h) {

    double sigma;

    if (D == 1) {
        sigma = 2./3.;
    }
    else if (D == 2) {
        sigma = 10./(7.*3.1415927);
    }
    else if (D == 3) {
        sigma = 1./3.1415927;
    }

    return 1./pow(h,D+2)*M_4(r/h)*sigma;
};

void SPH::calculate_rho(Particle *sources, int nsources, Particle *sinks, 
        int nsinks) {

    for (int i=0; i<nsinks; i++) {
        sinks[i].rho = 0.0;
        for (int j=0; j<nsources; j++) {
            double dr[3];
            VVV(dr, =sinks[i].r, -sources[j].r);
            
            double r;
            r = NORM(dr);

            sinks[i].rho += sources[j].m * kernel(r, sinks[i].h);
        }
    }
};

void SPH::calculate_accel(Particle *sources, int nsources, Particle *sinks, 
        int nsinks) {

    for (int i=0; i<nsinks; i++) {
        VS(sinks[i].a,=0);

        for (int j=0; j<nsources; i++) {
            double dr[3];
            VVV(dr, =sinks[i].r, -sources[j].r);
            
            double r;
            r = NORM(dr);

            double Fab = grad_kernel(r, sinks[i].h);

            double C = -sources[j].m*(sinks[i].P/(sinks[i].Omega*sinks[i].rho*
                        sinks[i].rho)*Fab + sources[j].P/(sources[j].Omega*
                        sources[j].rho*sources[j].rho)*Fab);

            VSV(sinks[i].a,+=C,*dr);

            VV(sinks[i].a,+=sinks[i].g);
        }
    }
};

void SPH::drift(double dt) {
    for (int i=0; i<S->np; i++)
        VV(P[i].r,+=dt*P[i].v);
};

void SPH::kick(double dt) {
    for (int i=0; i<S->np; i++)
        VV(P[i].v,+=dt*P[i].a);
};

void SPH::step() {
    for (int i=0; i<S->np; i++) {
        drift(S->dt/2.);

        //compute_accel(p, nparticles, p, nparticles);

        kick(S->dt/2.);

        update_timestep();

        //compute_accel(p, nparticles, p, nparticles);

        kick(S->dt);

        drift(S->dt);
    }
};

void SPH::update_timestep() {
};

double SPH::eos(double u, double rho) {
    double K=1.0;
    double gamma=5./3.;

    return K*pow(rho,gamma);
}
