/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains functions to assist problem initialization in init.c.

  The init_tools.c file contains helper functions used mainly in init.c. This is to keep init.c relatively clean, with only functions that were originally defined. 

  \author AYW 
  \date   2012-12-07 14:41 JST
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"

/* Global struct and arrays for normalization */
VarNorm vn;
double ini_cgs[USER_DEF_PARAMETERS];
double ini_code[USER_DEF_PARAMETERS];


/* Functions */

/* ************************************************************** */
void PrintInitData01(const double *prim) {
/*
 * Print some additional data during initialization
 *
 **************************************************************** */

    print1("\n");
    print1("> Conditions of ambient medium:\n");
    print1("\n");
    print1("  rho   %14g 1/cm3\n", prim[RHO] * vn.dens_norm / (MU_NORM * CONST_amu));
    print1("  te    %14g K\n", prim[PRS] / prim[RHO] * MU_NORM * KELVIN);
    print1("  cs    %14g km/s\n", sqrt(g_gamma * prim[PRS] / prim[RHO]) * vn.v_norm / 1.e5);
    print1("\n");
}


/* ************************************************************** */
void SetBaseNormalization() {
/*
 * Sets initializes VarNorm struct with derived normalizations
 * Gives cgs units upon multiplication form code units.
 *
 * Note that pot_norm is also the (mass) specific energy density [erg / g],
 * eint_norm is the energy [erg], and pres_norm is energy density [erg / cm^-3]
 *
 **************************************************************** */

    vn.l_norm = UNIT_LENGTH;
    vn.dens_norm = UNIT_DENSITY;
    vn.v_norm = UNIT_VELOCITY;
    vn.temp_norm = KELVIN;

    /* Derived normalizations */
    vn.t_norm = vn.l_norm / vn.v_norm;
    vn.area_norm = vn.l_norm * vn.l_norm;
    vn.pres_norm = vn.dens_norm * vn.v_norm * vn.v_norm;
    vn.power_norm = vn.pres_norm * vn.v_norm * vn.area_norm;
    vn.eflux_norm = vn.pres_norm * vn.v_norm;
    vn.eint_norm = vn.pres_norm * vn.l_norm * vn.l_norm * vn.l_norm;
    vn.mdot_norm = vn.dens_norm * pow(vn.l_norm, 3) / vn.t_norm;
    vn.newton_norm = 1. / (vn.t_norm * vn.t_norm * vn.dens_norm);
    vn.pot_norm = vn.v_norm * vn.v_norm;
    vn.acc_norm = vn.v_norm / vn.t_norm;
    vn.n_norm = 1. / (vn.l_norm * vn.l_norm * vn.l_norm);
    vn.m_norm = vn.dens_norm * vn.l_norm * vn.l_norm * vn.l_norm;

    print1("> Base normalization initialized.\n\n");

    return;
}


/* ************************************************************** */
void SetIniNormalization() {
/*
 * Sets noramlizations for ini input variables.
 * ini_cgs converts ini parameters to cgs upon multiplication,
 * whereas ini_code converts ini parameters to code units.
 *
 **************************************************************** */

    double year, degrad;
    year = CONST_ly / CONST_c;
    degrad = CONST_PI / 180.;

    ini_cgs[PAR_MACH] = 1.;
    ini_cgs[PAR_DENS] = 1.;
    ini_cgs[PAR_DRATIO] = 1.;

    ini_code[PAR_MACH] = 1.;
    ini_code[PAR_DENS] = 1.;
    ini_code[PAR_DRATIO] = 1.;


    print1("> Ini parameter normalization array initialized.\n\n");

    return;

}


/* ************************************************ */
void PrintGridStruct(Grid *grid, int show_for_rank, int k, int j, int i) {
    /*!
     * grid     array of grid structures
     *
     * The function prints out grid structure members and
     * is useful for parallel debugging.
     *
     ************************************************** */

    if (prank == show_for_rank) {
        printf("\n");
        printf("Printing members of GRID structure on rank %d\n", prank);
        printf("                        IDIR             JDIR              KDIR\n");
        printf("grid[].xi          = %16e %16e %16e\n", grid[IDIR].xi,          grid[JDIR].xi,          grid[KDIR].xi);
        printf("grid[].xf          = %16e %16e %16e\n", grid[IDIR].xf,          grid[JDIR].xf,          grid[KDIR].xf);
        printf("grid[].x[]         = %16e %16e %16e\n", grid[IDIR].x[i],        grid[JDIR].x[j],        grid[KDIR].x[k]);
        printf("grid[].x_glob[]    = %16e %16e %16e\n", grid[IDIR].x_glob[i],   grid[JDIR].x_glob[j],   grid[KDIR].x_glob[k]);
        printf("grid[].xr[]        = %16e %16e %16e\n", grid[IDIR].xr[i],       grid[JDIR].xr[j],       grid[KDIR].xr[k]);
        printf("grid[].xr_glob[]   = %16e %16e %16e\n", grid[IDIR].xr_glob[i],  grid[JDIR].xr_glob[j],  grid[KDIR].xr_glob[k]);
        printf("grid[].xl[]        = %16e %16e %16e\n", grid[IDIR].xl[i],       grid[JDIR].xl[j],       grid[KDIR].xl[k]);
        printf("grid[].xl_glob[]   = %16e %16e %16e\n", grid[IDIR].xl_glob[i],  grid[JDIR].xl_glob[j],  grid[KDIR].xl_glob[k]);
        printf("grid[].dx[]        = %16e %16e %16e\n", grid[IDIR].dx[i],       grid[JDIR].dx[j],       grid[KDIR].dx[k]);
        printf("grid[].dx_glob[]   = %16e %16e %16e\n", grid[IDIR].dx_glob[i],  grid[JDIR].dx_glob[j],  grid[KDIR].dx_glob[k]);
        printf("grid[].xgc[]       = %16e %16e %16e\n", grid[IDIR].xgc[i],      grid[JDIR].xgc[j],      grid[KDIR].xgc[k]);
        printf("grid[].dV[]        = %16e %16e %16e\n", grid[IDIR].dV[i],       grid[JDIR].dV[j],       grid[KDIR].dV[k]);
        printf("grid[].A[]         = %16e %16e %16e\n", grid[IDIR].A[i],        grid[JDIR].A[j],        grid[KDIR].A[k]);
        printf("grid[].r_1[]       = %16e %16e %16e\n", grid[IDIR].r_1[i],      grid[JDIR].r_1[j],      grid[KDIR].r_1[k]);
        printf("grid[].ct[]        = %16e %16e %16e\n", grid[IDIR].ct[i],       grid[JDIR].ct[j],       grid[KDIR].ct[k]);
        printf("grid[].inv_dx[]    = %16e %16e %16e\n", grid[IDIR].inv_dx[i],   grid[JDIR].inv_dx[j],   grid[KDIR].inv_dx[k]);
        printf("grid[].inv_dxi[]   = %16e %16e %16e\n", grid[IDIR].inv_dxi[i],  grid[JDIR].inv_dxi[j],  grid[KDIR].inv_dxi[k]);
        printf("grid[].dl_min      = %16e %16e %16e\n", grid[IDIR].dl_min,      grid[JDIR].dl_min,      grid[KDIR].dl_min);
        printf("grid[].np_tot_glob = %16d %16d %16d\n", grid[IDIR].np_tot_glob, grid[JDIR].np_tot_glob, grid[KDIR].np_tot_glob);
        printf("grid[].np_int_glob = %16d %16d %16d\n", grid[IDIR].np_int_glob, grid[JDIR].np_int_glob, grid[KDIR].np_int_glob);
        printf("grid[].np_tot      = %16d %16d %16d\n", grid[IDIR].np_tot,      grid[JDIR].np_tot,      grid[KDIR].np_tot);
        printf("grid[].np_int      = %16d %16d %16d\n", grid[IDIR].np_int,      grid[JDIR].np_int,      grid[KDIR].np_int);
        printf("grid[].nghost      = %16d %16d %16d\n", grid[IDIR].nghost,      grid[JDIR].nghost,      grid[KDIR].nghost);
        printf("grid[].lbound      = %16d %16d %16d\n", grid[IDIR].lbound,      grid[JDIR].lbound,      grid[KDIR].lbound);
        printf("grid[].rbound      = %16d %16d %16d\n", grid[IDIR].rbound,      grid[JDIR].rbound,      grid[KDIR].rbound);
        printf("grid[].gbeg        = %16d %16d %16d\n", grid[IDIR].gbeg,        grid[JDIR].gbeg,        grid[KDIR].gbeg);
        printf("grid[].gend        = %16d %16d %16d\n", grid[IDIR].gend,        grid[JDIR].gend,        grid[KDIR].gend);
        printf("grid[].beg         = %16d %16d %16d\n", grid[IDIR].beg,         grid[JDIR].beg,         grid[KDIR].beg);
        printf("grid[].end         = %16d %16d %16d\n", grid[IDIR].end,         grid[JDIR].end,         grid[KDIR].end);
        printf("grid[].lbeg        = %16d %16d %16d\n", grid[IDIR].lbeg,        grid[JDIR].lbeg,        grid[KDIR].lbeg);
        printf("grid[].lend        = %16d %16d %16d\n", grid[IDIR].lend,        grid[JDIR].lend,        grid[KDIR].lend);
        printf("grid[].uniform     = %16d %16d %16d\n", grid[IDIR].uniform,     grid[JDIR].uniform,     grid[KDIR].uniform);
        printf("grid[].nproc       = %16d %16d %16d\n", grid[IDIR].nproc,       grid[JDIR].nproc,       grid[KDIR].nproc);
        printf("grid[].rank_coord  = %16d %16d %16d\n", grid[IDIR].rank_coord,  grid[JDIR].rank_coord,  grid[KDIR].rank_coord);
        printf("grid[].level       = %16d %16d %16d\n", grid[IDIR].level,       grid[JDIR].level,       grid[KDIR].level);

    }

}

/* ************************************************ */
void PrintBaseNormalizations() {
    /*!
     * grid     array of grid structures
     *
     * The function prints out grid structure members and
     * is useful for parallel debugging.
     *
     ************************************************** */

    print1("vn.l_norm      = %16e \n", vn.l_norm );
    print1("vn.dens_norm   = %16e \n", vn.dens_norm );
    print1("vn.v_norm      = %16e \n", vn.v_norm );
    print1("vn.temp_norm   = %16e \n", vn.temp_norm );
    print1("vn.t_norm      = %16e \n", vn.t_norm );
    print1("vn.area_norm   = %16e \n", vn.area_norm );
    print1("vn.pres_norm   = %16e \n", vn.pres_norm );
    print1("vn.power_norm  = %16e \n", vn.power_norm );
    print1("vn.eflux_norm  = %16e \n", vn.eflux_norm );
    print1("vn.eint_norm   = %16e \n", vn.eint_norm );
    print1("vn.mdot_norm   = %16e \n", vn.mdot_norm );
    print1("vn.newton_norm = %16e \n", vn.newton_norm );
    print1("vn.pot_norm    = %16e \n", vn.pot_norm );
    print1("vn.acc_norm    = %16e \n", vn.acc_norm );
    print1("vn.n_norm      = %16e \n", vn.n_norm );
    print1("vn.m_norm      = %16e \n", vn.m_norm );
    print1("\n");
    print1("\n");
    print1("\n");

}
