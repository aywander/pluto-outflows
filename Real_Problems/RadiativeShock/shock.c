//
// Created by Alexander Y. Wagner on 10/30/16.
//

#include "shock.h"
#include "pluto.h"

double mach_from_te1(const double spd, const double te, const double mu) {

    return sqrt( -(4 * (g_gamma - 1) * spd * spd ) /
                 ((1 + (g_gamma - 6) * g_gamma) * spd * spd +
                  (1 + g_gamma) * spd * sqrt(-8 * (g_gamma - 1) * te * g_gamma / mu +
                                                     pow(1 + g_gamma, 2) * spd * spd)));
}

double te_jump(const double mach) {

    return (2 * g_gamma * mach * mach - (g_gamma - 1)) * ((g_gamma - 1) * mach * mach + 2) / \
               pow((g_gamma + 1) * mach, 2);

}

double rho_jump(const double mach) {

    return (g_gamma + 1) * mach * mach / ((g_gamma - 1) * mach * mach + 2);

}

double rho_jump_cdl(const double mach, const double t0, const double td) {

    double tr = t0 / td;

    double p = (g_gamma * mach * mach + 1.) * tr;
    double q = g_gamma * mach * mach * tr;

    return p / 2. + sqrt(p * p / 4. - q);
    // return g_gamma * mach * mach; // This is only in the isothermal and strong shock limit

}

double mach(const double u, const double te, const double mu) {
    return u / sqrt(g_gamma * te / mu);
}

ShockState sh;
