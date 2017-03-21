#ifndef PLUTO_USR_H
#define PLUTO_USR_H

/* Shock initialization methods -- SHOCK_INIT_MODE */
#define SI_IMPULSIVE    0
#define SI_WALL         1

/* Shock initial condition calculations -- SHOCK_COND_MODE */
#define SC_UPSTREAM     0
#define SC_SHOCK_TE     1

/* MU_CALC values. Methods of calculating the mean molecular mass. */
#define MU_CONST        0
#define MU_TABLE        1
#define MU_ANALYTIC     2
#define MU_FRACTIONS    3
#define MU_FNAME "mutable.dat"


/* Include user definitions and user macros here (last) */
#include "definitions_usr.h"
#include "macros_usr.h"

#endif
