#include <cmath>
#include <stdio.h>
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

template<int D>
struct Params {
    double rmax[D];
    double time;
    double dt;
    double nsteps;
    double eta;
    double alpha_min;
    double beta;
    double epsilon;
    double sigma;
    double kernel_const;
    int np;
    int Lboundary[D];
    int Rboundary[D];
};

template<int D>
struct SPH{
    Particle<D> *P;
    Params<D> *S;

    // Math helper functions.
    double kernel(double r, double h);
    double grad_kernel(double r, double h);
    double kernel_hderiv(double r, double h);
    double eos(double u, double rho);
    double sound_speed(double P, double rho);
    double del_dot_v(Particle<D> *sink, Particle<D> *sources, int nsources);
    double laplace_u(Particle<D> *sink, Particle<D> *sources, int nsources);
    Vector<double, D> delta_r(Particle<D> *a, Particle<D> *b);

    // Functions to calculate elements of the Particle structure.
    void calculate_P();
    void calculate_cs();
    void calculate_h(Particle<D> *sinks, int nsinks, Particle<D> *sources, 
            int nsources);
    void calculate_Omega(Particle<D> *sinks, int nsinks, Particle<D> *sources,
            int nsources);
    void calculate_rho(Particle<D> *sinks, int nsinks, Particle<D> *sources, 
            int nsources);
    void calculate_accel(Particle<D> *sinks, int nsinks, Particle<D> *sources,
            int nsources);
    void calculate_uvel(Particle<D> *sinks, int nsinks, Particle<D> *sources,
            int nsources);
    void calculate_alphavel(Particle<D> *sinks, int nsinks, 
            Particle<D> *sources, int nsources);
    void calculate_alphauvel(Particle<D> *sinks, int nsinks, 
            Particle<D> *sources, int nsources);
    void calculate_diss(Particle<D> *sinks, int nsinks, Particle<D> *sources,
            int nsources);
    void calculate_udiss(Particle<D> *sinks, int nsinks, Particle<D> *sources,
            int nsources);

    // Functions to move particles.
    void drift(double dt);
    void drift_u(double dt);
    void drift_alpha(double dt);
    void drift_alphau(double dt);
    void kick(double dt);
    void apply_boundary_conditions();
    void update_timestep();
    void initialize();
    void first_step();
    void step();
};

// Various useful math routines.

template<int D>
double SPH<D>::kernel(double r, double h) {

    double q = r/h;

    if ((0 <= q) && (q < 1))
        return (1./pow(h,D))*
            (1./4.*(2.-q)*(2.-q)*(2.-q)-(1.-q)*(1.-q)*(1.-q))*S->kernel_const;
    else if ((1 <= q) && (q <= 2))
        return (1./pow(h,D))*(1./4.*(2-q)*(2-q)*(2-q))*S->kernel_const;
    else
        return 0.0;
};

template<int D>
double SPH<D>::grad_kernel(double r, double h) {

    double q = r/h;

    if ((0 <= q) && (q < 1))
        return (1./pow(h,D+2))*(-3.+9./4.*q)*S->kernel_const;
    else if ((1 <= q) && (q <= 2))
        return (1./pow(h,D+2))*(-3./4.*(2.-q)*(2.-q)/q)*S->kernel_const;
    else
        return 0.0;
};

template<int D>
double SPH<D>::kernel_hderiv(double r, double h) {

    double q = r/h;

    if ((0 <= q) && (q < 1))
        return (1./pow(h,D+1))*(-D+(3.+3./2.*D)*q*q-(9./4.+3./4.*D)*q*q*q)*
            S->kernel_const;
    else if ((1 <= q) && (q <= 2))
        return (1./pow(h,D+1))*(-2.*D+(3.+3.*D)*q-(3.+3./2.*D)*q*q+
            (3./4.+1./4.*D)*q*q*q)*S->kernel_const;
    else
        return 0.0;
}

template<int D>
double SPH<D>::del_dot_v(Particle<D> *sink, Particle<D> *sources, int nsources) {
    double div_v = 0.0;

    for (int j=0; j<nsources; j++) {
        Vector<double, D> dv = sink->v - sources[j].v;

        Vector<double, D> dr = delta_r(sink,&sources[j]);

        double r = dr.norm();
            
        double dv_dot_dr = dv*dr;

        double Fab = grad_kernel(r,sink->h);

        div_v += -sources[j].m/sources[j].rho*dv_dot_dr*Fab;
    }

    return div_v;
}

template<int D>
double SPH<D>::laplace_u(Particle<D> *sink, Particle<D> *sources, 
        int nsources) {
    double du = 0.0;

    for (int j=0; j<nsources; j++) {
        Vector<double, D> dr = delta_r(sink,&sources[j]);

        double r = dr.norm();

        double Fab = grad_kernel(r,sink->h);

        du += sources[j].m/sources[j].rho*(sink->u-sources[j].u)*Fab;
    }

    return 2*du;
}

template<int D>
Vector<double, D> SPH<D>::delta_r(Particle<D> *a, Particle<D> *b) {
    Vector<double,D> dr = a->r - b->r;

    if ((S->Lboundary[0] == 1) && (S->Rboundary[0] == 1)) {
        for (int k=0; k<D; k++) {
            if (dr[k] > 0.5)
                dr[k] -= 1.0;
            else if (dr[k] < -0.5*S->rmax[k])
                dr[k] += 1.0;
        }
    }
    
    return dr;
}

template<int D>
double SPH<D>::eos(double u, double rho) {
    double gamma=7./5.;

    return (gamma-1.0)*rho*u;
}

template<int D>
double SPH<D>::sound_speed(double P, double rho) {
    double gamma=7./5.;

    return sqrt(gamma*P/rho);
}

// Routines to calculate the various parameters of the Particle class.

template<int D>
void SPH<D>::calculate_P() {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].P = eos(P[i].u, P[i].rho);
};

template<int D>
void SPH<D>::calculate_cs() {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].cs = sound_speed(P[i].P, P[i].rho);
}

template<int D>
void SPH<D>::calculate_Omega(Particle<D> *sources, int nsources, 
        Particle<D> *sinks, int nsinks) {
#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        sinks[i].Omega = 1.0;

        double C = -sinks[i].h/(D*sinks[i].rho);

        for (int j=0; j<nsources; j++) {
            Vector<double, D> dr = delta_r(&sinks[i],&sources[j]);
            double r = dr.norm();

            sinks[i].Omega -= C*sources[j].m*kernel_hderiv(r,sinks[i].h);
        }
    }
}

template<int D>
void SPH<D>::calculate_h(Particle<D> *sources, int nsources, Particle<D> *sinks,
        int nsinks) {

    bool not_converged = true;

    while (not_converged) {
        calculate_rho(sources, nsources, sinks, nsinks);
        calculate_Omega(sources, nsources, sinks, nsinks);

        int converge_test = 0;

        for (int i=0; i<nsinks; i++) {
            double rho_ha = pow(S->eta,D)*sinks[i].m/pow(sinks[i].h,D);
            double fprime = -D*sinks[i].rho*sinks[i].Omega/sinks[i].h;

            double hold = sinks[i].h;
            sinks[i].h -= (rho_ha-sinks[i].rho)/fprime;

            double diff = fabs(hold-sinks[i].h)/hold;

            if (diff < 0.001)
                converge_test++;
        }

        if (converge_test == nsinks)
            not_converged = false;
    }
}

template<int D>
void SPH<D>::calculate_rho(Particle<D> *sources, int nsources, 
        Particle<D> *sinks, int nsinks) {

#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        sinks[i].rho = 0.0;
        for (int j=0; j<nsources; j++) {
            Vector<double, D> dr = delta_r(&sinks[i],&sources[j]);

            double r = dr.norm();

            sinks[i].rho += sources[j].m * kernel(r, sinks[i].h);
        }
    }
};

template<int D>
void SPH<D>::calculate_accel(Particle<D> *sources, int nsources, 
        Particle<D> *sinks, int nsinks) {

#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        sinks[i].a = 0;

        for (int j=0; j<nsources; j++) {
            Vector<double, D> dr = delta_r(&sinks[i],&sources[j]);
            
            double r = dr.norm();

            double Fabha = grad_kernel(r, sinks[i].h);
            double Fabhb = grad_kernel(r, sources[j].h);

            double C = -sources[j].m*(sinks[i].P/(sinks[i].Omega*sinks[i].rho*
                        sinks[i].rho)*Fabha + sources[j].P/(sources[j].Omega*
                        sources[j].rho*sources[j].rho)*Fabhb);

            sinks[i].a += C*dr;
        }
    }
};

template<int D>
void SPH<D>::calculate_uvel(Particle<D> *sinks, int nsinks, 
        Particle<D> *sources, int nsources) {
    
#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        sinks[i].vu = 0.0;

        double C = sinks[i].P / (sinks[i].Omega * sinks[i].rho * sinks[i].rho);

        for (int j=0; j<nsources; j++) {
            Vector<double, D> dr = delta_r(&sinks[i],&sources[j]);
            
            double r = dr.norm();

            double Fab = grad_kernel(r, sinks[i].h);

            Vector<double, D> dv = sinks[i].v - sources[j].v;

            double dv_dot_dr = dv*dr;

            sinks[i].vu += C*sources[j].m*dv_dot_dr*Fab;
        }
    }
}

template<int D>
void SPH<D>::calculate_alphavel(Particle<D> *sinks, int nsinks, 
        Particle<D> *sources, int nsources) {

#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        double dv = del_dot_v(&sinks[i], sources, nsources);

        double tau = sinks[i].h/(S->sigma*sinks[i].cs);

        double s = 0.0;
        if (-dv > 0.0)
            s=-dv;

        sinks[i].valpha = s - (sinks[i].alpha-S->alpha_min)/tau;
    }
}

template<int D>
void SPH<D>::calculate_alphauvel(Particle<D> *sinks, int nsinks, 
        Particle<D> *sources, int nsources) {

#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        double du = laplace_u(&sinks[i], sources, nsources);

        double tau = sinks[i].h/(S->sigma*sinks[i].cs);

        double s = fabs(0.1*sinks[i].h*du);

        sinks[i].valphau = s - (sinks[i].alphau-S->alpha_min)/tau;
    }
}

template<int D>
void SPH<D>::calculate_diss(Particle<D> *sinks, int nsinks, 
        Particle<D> *sources, int nsources) {
#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        sinks[i].adiss = 0;

        for (int j=0; j<nsources; j++) {
            Vector<double, D> dv = sinks[i].v - sources[j].v;

            Vector<double, D> dr = delta_r(&sinks[i],&sources[j]);
            
            double rhoab = 0.5*(sinks[i].rho+sources[j].rho);

            double r = dr.norm();
            double Fab = 0.5*(grad_kernel(r,sinks[i].h) + 
                    grad_kernel(r,sources[j].h));

            double dv_dot_dr = dv*dr;

            double vsig = 0.0;
            if (dv_dot_dr <= 0) {
                if (r == 0.0)
                    vsig = 0.5*(sinks[i].cs + sources[j].cs);
                else
                    vsig = 0.5*(sinks[i].cs + sources[j].cs - 
                            S->beta*dv_dot_dr/r);
            }

            double C = 0.0;
            if (r != 0.0)
                C = -sources[j].m*(-sinks[i].alpha*vsig*dv_dot_dr/
                        (rhoab*r))*Fab;

            sinks[i].adiss += C*dr;
        }
    }
}

template<int D>
void SPH<D>::calculate_udiss(Particle<D> *sinks, int nsinks, 
        Particle<D> *sources, int nsources) {
#pragma omp parallel for
    for (int i=0; i<nsinks; i++) {
        sinks[i].vudiss = 0.0;
        for (int j=0; j<nsources; j++) {
            Vector<double, D> dr = delta_r(&sinks[i],&sources[j]);
            double r = dr.norm();

            Vector<double, D> dv = sinks[i].v - sources[j].v;

            double rhoab = 0.5*(sinks[i].rho+sources[j].rho);

            double vusig = sqrt(fabs(sinks[i].P-sources[j].P)/rhoab);
            
            double uab = sinks[i].u - sources[j].u;

            double Fab = 0.5*(grad_kernel(r,sinks[i].h) + 
                    grad_kernel(r,sources[j].h));

            double dv_dot_dr = dv*dr;

            double vsig = 0.0;
            if (dv_dot_dr <= 0.0) {
                if (r == 0.0)
                    vsig = 0.5*(sinks[i].cs + sources[j].cs);
                else
                    vsig = 0.5*(sinks[i].cs + sources[j].cs - 
                            S->beta*dv_dot_dr/r);
            }

            double C1 = 0.0;
            if (r != 0.0)
                C1 = -0.5*sinks[i].alpha*vsig*dv_dot_dr*
                    dv_dot_dr/(r*r);
            double C2 = sinks[i].alphau*vusig*uab;

            sinks[i].vudiss += sources[j].m/rhoab*(C1+C2)*r*Fab;
        }
    }
}

// Routines to move the particles.

template<int D>
void SPH<D>::drift(double dt) {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].r += dt * P[i].v;
    apply_boundary_conditions();
};

template<int D>
void SPH<D>::drift_u(double dt) {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].u += dt*(P[i].vu+P[i].vudiss);
}

template<int D>
void SPH<D>::drift_alpha(double dt) {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].alpha += dt*P[i].valpha;
}

template<int D>
void SPH<D>::drift_alphau(double dt) {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].alphau += dt*P[i].valphau;
}

template<int D>
void SPH<D>::kick(double dt) {
#pragma omp parallel for
    for (int i=0; i<S->np; i++)
        P[i].v += dt * P[i].a + dt * P[i].adiss;
};

template<int D>
void SPH<D>::update_timestep() {
    double test;
    for (int i=0; i<S->np; i++) {
        double dv = del_dot_v(&P[i], P, S->np);

        double v = P[i].h*fabs(dv)+P[i].cs+1.2*P[i].alpha*P[i].cs;
        if (dv < 0)
            v += S->beta*P[i].h*fabs(dv);

        double temp_dt = 0.25*P[i].h/v;

        if (i == 0)
            S->dt = temp_dt;
        else if (temp_dt < S->dt)
            S->dt = temp_dt;
    }
};

template<int D>
void SPH<D>::apply_boundary_conditions() {
    for (int i=0; i<S->np; i++) {
        //
        // 0 = No boundary conditions, 1 = Periodic, 2 = Reflecting
        //
        for (int j=0; j<D; j++) {
            if (S->Lboundary[j] == 1) {
                P[i].r[j] += S->rmax[j];
                P[i].r[j] = fmod(P[i].r[j],S->rmax[j]);
            }
            else if ((S->Lboundary[j] == 2) && (P[i].r[j] <= 0.0)) {
                P[i].r[j] *= -1;
                P[i].v[j] *= -1;
                P[i].a[j] *= -1;
            }
        }

        for (int j=0; j<D; j++) {
            if (S->Lboundary[j] == 1) {
                P[i].r[j] += S->rmax[j];
                P[i].r[j] = fmod(P[i].r[j],S->rmax[j]);
            }
            else if ((S->Lboundary[j] == 2) && (P[i].r[j] > 1.0)) {
                P[i].r[j] = 2-P[i].r[j];
                P[i].v[j] *= -1;
                P[i].a[j] *= -1;
            }
        }
    }
}

template<int D>
void SPH<D>::initialize() {
    calculate_h(P, S->np, P, S->np);
    calculate_rho(P, S->np, P, S->np);
    calculate_Omega(P, S->np, P, S->np);
    calculate_P();
    calculate_cs();
}

template<int D>
void SPH<D>::step() {
    update_timestep();

    drift(S->dt/2.);
    drift_u(S->dt/2.);
    drift_alpha(S->dt/2.);
    drift_alphau(S->dt/2.);

    calculate_h(P, S->np, P, S->np);
    calculate_rho(P, S->np, P, S->np);
    calculate_Omega(P, S->np, P, S->np);
    calculate_P();
    calculate_cs();

    calculate_accel(P, S->np, P, S->np);
    calculate_diss(P, S->np, P, S->np);

    kick(S->dt);

    calculate_uvel(P, S->np, P, S->np);
    calculate_udiss(P, S->np, P, S->np);
    calculate_alphavel(P, S->np, P, S->np);
    calculate_alphauvel(P, S->np, P, S->np);

    drift(S->dt/2.);
    drift_u(S->dt/2.);
    drift_alpha(S->dt/2.);
    drift_alphau(S->dt/2.);

    S->time += S->dt;
};

// wrappers for Python

extern "C" {
    void sph_new_1D(SPH<1> *me, Params<1> *S, Particle<1> *P) {
        me->S = S;
        me->P = P;
    }
    void sph_new_2D(SPH<2> *me, Params<2> *S, Particle<2> *P) {
        me->S = S;
        me->P = P;
    }
    void sph_new_3D(SPH<3> *me, Params<3> *S, Particle<3> *P) {
        me->S = S;
        me->P = P;
    }

    void sph_step_1D(SPH<1> *me) {
        me->step();
    }
    void sph_step_2D(SPH<2> *me) {
        me->step();
    }
    void sph_step_3D(SPH<3> *me) {
        me->step();
    }

    void sph_initialize_1D(SPH<1> *me) {
        me->initialize();
    }
    void sph_initialize_2D(SPH<2> *me) {
        me->initialize();
    }
    void sph_initialize_3D(SPH<3> *me) {
        me->initialize();
    }
}
