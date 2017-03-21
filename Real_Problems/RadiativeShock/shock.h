//
// Created by Alexander Y. Wagner on 10/30/16.
//

#ifndef PLUTO_SHOCK_H
#define PLUTO_SHOCK_H

/* Shock state variables.
 *
 * labels:
 *   0 : distant upstream
 *   1 : post-(sub)shock
 *   p : pre-(sub)shock (precursor)
 *   d : distant downstream
 *
 * vars:
 *   xshock : location of shock (only in 1D)
 *
 * Note, assume we are in frame of the shock => Shock velocity = upstream velocity = u0.
 * If there is no precursor, <>_p = <>_0
 * If shock is non-radiative, <>_d = <>_1
 * */

typedef struct {
    double rho;
    double u;
    double te;
    double mu;
    double bprp;
    double bpar;
} FluidState;

typedef struct {

    FluidState s0;
    FluidState s1;
    FluidState sp;
    FluidState sd;
    double xshock;

} ShockState;

extern ShockState sh;


double mach_from_te1(const double spd, const double te, const double mu);

double te_jump(const double mach);

double rho_jump(const double mach);

double rho_jump_cdl(const double mach, const double t0, const double td);

double mach(const double u, const double te, const double mu);

#endif //PLUTO_SHOCK_H
