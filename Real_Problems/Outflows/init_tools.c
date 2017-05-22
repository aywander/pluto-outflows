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
void PrintInitData01(const double *out_primitives,
                     const double *halo_primitives) {
/*
 * Print some additional data during initialization
 *
 **************************************************************** */


    print1("\n");
    print1("> Conditions at (0, 0, 0):\n");
    print1("\n");
    print1("        %14s  %14s  %14s\n",
           "Nozzle primitives", "Hot halo prim", "ratio");
    print1("  rho   %14g  %14g  %14g\n",
           out_primitives[RHO], halo_primitives[RHO],
           out_primitives[RHO] / halo_primitives[RHO]);
    print1("  pr    %14g  %14g  %14g\n",
           out_primitives[PRS], halo_primitives[PRS],
           out_primitives[PRS] / halo_primitives[PRS]);
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

    ini_cgs[PAR_OPOW] = 1.;
    ini_cgs[PAR_OSPD] = NOZZLE_SELECT(1., CONST_c);
    ini_cgs[PAR_OMDT] = NOZZLE_SELECT(1., CONST_Msun / year);
    ini_cgs[PAR_OANG] = degrad;
    ini_cgs[PAR_ORAD] = vn.l_norm;
    ini_cgs[PAR_ODBH] = vn.l_norm;
    ini_cgs[PAR_ODIR] = degrad;
    ini_cgs[PAR_OOMG] = 1. / (1.e6 * year);
    ini_cgs[PAR_OPHI] = degrad;
    ini_cgs[PAR_OSPH] = vn.l_norm;
    ini_cgs[PAR_ARAD] = vn.l_norm;
    ini_cgs[PAR_AMBH] = CONST_Msun;
    ini_cgs[PAR_AEFF] = 1.;
    ini_cgs[PAR_ASNK] = vn.l_norm;
    ini_cgs[PAR_HRHO] = vn.dens_norm;
    ini_cgs[PAR_HTMP] = 1;
    ini_cgs[PAR_HVX1] = 1.e5;
    ini_cgs[PAR_HVX2] = 1.e5;
    ini_cgs[PAR_HVX3] = 1.e5;
    ini_cgs[PAR_HVRD] = 1.e5;
    ini_cgs[PAR_HRAD] = vn.l_norm;
//    ini_cgs[PAR_MGAL] = CONST_Msun;
    ini_cgs[PAR_WRHO] = vn.dens_norm;
    ini_cgs[PAR_WTRB] = 1.e5;
    ini_cgs[PAR_WRAD] = vn.l_norm;
    ini_cgs[PAR_WROT] = 1.;
    ini_cgs[PAR_WX1L] = vn.l_norm;
    ini_cgs[PAR_WX1H] = vn.l_norm;
    ini_cgs[PAR_WX2L] = vn.l_norm;
    ini_cgs[PAR_WX2H] = vn.l_norm;
    ini_cgs[PAR_WX3L] = vn.l_norm;
    ini_cgs[PAR_WX3H] = vn.l_norm;
    ini_cgs[PAR_WVRD] = 1.e5;
    ini_cgs[PAR_WVPL] = 1.e5;
    ini_cgs[PAR_WVPP] = 1.e5;
    ini_cgs[PAR_WVAN] = 1.e5;
    ini_cgs[PAR_SGAV] = 1.e5;
    ini_cgs[PAR_NCLD] = 1;
//    ini_cgs[PAR_SPOW] = 1;
//    ini_cgs[PAR_SDUR] = 1.e6 * year;
//    ini_cgs[PAR_SRAD] = vn.l_norm;
//    ini_cgs[PAR_SHGT] = vn.l_norm;
    ini_cgs[PAR_LOMX] = 1;
    ini_cgs[PAR_LCMX] = 1;

    ini_code[PAR_OPOW] = ini_cgs[PAR_OPOW] / vn.power_norm;
    ini_code[PAR_OSPD] = ini_cgs[PAR_OSPD] / NOZZLE_SELECT(1., vn.v_norm);
    ini_code[PAR_OMDT] = ini_cgs[PAR_OMDT] / NOZZLE_SELECT(1., vn.mdot_norm);
    ini_code[PAR_OANG] = ini_cgs[PAR_OANG];
    ini_code[PAR_ORAD] = ini_cgs[PAR_ORAD] / vn.l_norm;
    ini_code[PAR_ODBH] = ini_cgs[PAR_ODBH] / vn.l_norm;
    ini_code[PAR_ODIR] = ini_cgs[PAR_ODIR];
    ini_code[PAR_OOMG] = ini_cgs[PAR_OOMG] * vn.t_norm;
    ini_code[PAR_OPHI] = ini_cgs[PAR_OPHI];
    ini_code[PAR_OSPH] = ini_cgs[PAR_OSPH] / vn.l_norm;
    ini_code[PAR_ARAD] = ini_cgs[PAR_ARAD] / vn.l_norm;
    ini_code[PAR_AMBH] = ini_cgs[PAR_AMBH] / vn.m_norm;
    ini_code[PAR_AEFF] = ini_cgs[PAR_AEFF];
    ini_code[PAR_ASNK] = ini_cgs[PAR_ASNK] / vn.l_norm;
    ini_code[PAR_HRHO] = ini_cgs[PAR_HRHO] / vn.dens_norm;
    ini_code[PAR_HTMP] = ini_cgs[PAR_HTMP] / vn.temp_norm;
    ini_code[PAR_HVX1] = ini_cgs[PAR_HVX1] / vn.v_norm;
    ini_code[PAR_HVX2] = ini_cgs[PAR_HVX2] / vn.v_norm;
    ini_code[PAR_HVX3] = ini_cgs[PAR_HVX3] / vn.v_norm;
    ini_code[PAR_HVRD] = ini_cgs[PAR_HVRD] / vn.v_norm;
    ini_code[PAR_HRAD] = ini_cgs[PAR_HRAD] / vn.l_norm;
//    ini_code[PAR_MGAL] = ini_cgs[PAR_MGAL] / vn.m_norm;
    ini_code[PAR_WRHO] = ini_cgs[PAR_WRHO] / vn.dens_norm;
    ini_code[PAR_WTRB] = ini_cgs[PAR_WTRB] / vn.v_norm;
    ini_code[PAR_WRAD] = ini_cgs[PAR_WRAD] / vn.l_norm;
    ini_code[PAR_WROT] = ini_cgs[PAR_WROT];
    ini_code[PAR_WX1L] = ini_cgs[PAR_WX1L] / vn.l_norm;
    ini_code[PAR_WX1H] = ini_cgs[PAR_WX1H] / vn.l_norm;
    ini_code[PAR_WX2L] = ini_cgs[PAR_WX2L] / vn.l_norm;
    ini_code[PAR_WX2H] = ini_cgs[PAR_WX2H] / vn.l_norm;
    ini_code[PAR_WX3L] = ini_cgs[PAR_WX3L] / vn.l_norm;
    ini_code[PAR_WX3H] = ini_cgs[PAR_WX3H] / vn.l_norm;
    ini_code[PAR_WVRD] = ini_cgs[PAR_WVRD] / vn.v_norm;
    ini_code[PAR_WVPL] = ini_cgs[PAR_WVPL] / vn.v_norm;
    ini_code[PAR_WVPP] = ini_cgs[PAR_WVPP] / vn.v_norm;
    ini_code[PAR_WVAN] = ini_cgs[PAR_WVAN] / vn.v_norm;
    ini_code[PAR_SGAV] = ini_cgs[PAR_SGAV] / vn.v_norm;
    ini_code[PAR_NCLD] = ini_cgs[PAR_NCLD];
//    ini_code[PAR_SPOW] = ini_cgs[PAR_SPOW] / vn.power_norm;
//    ini_code[PAR_SDUR] = ini_cgs[PAR_SDUR] / vn.t_norm;
//    ini_code[PAR_SRAD] = ini_cgs[PAR_SRAD] / vn.l_norm;
//    ini_code[PAR_SHGT] = ini_cgs[PAR_SHGT] / vn.l_norm;
    ini_code[PAR_LOMX] = ini_cgs[PAR_LOMX];
    ini_code[PAR_LCMX] = ini_cgs[PAR_LCMX];

//----DM: 7Aug15, Turned off PAR_LEV1 as not enough parameter slots---//
    /* AYW 2016-02-29 NOTE: increased number of parameters */
//  ini_cgs[PAR_LEV1] = 1;
//  ini_code[PAR_LEV1] = ini_cgs[PAR_LEV1];


    print1("> Ini parameter normalization array initialized.\n\n");

    return;

}


/* ************************************************ */
double Speed2Lorentz(const double vel)
/*!
 * Return Lorentz factor from Lorentz factor
 *
 ************************************************** */
{
    double vsqr = 0, clight2;

    clight2 = pow(CONST_c / UNIT_VELOCITY, 2);

    vsqr = vel * vel;
    return 1.0 / sqrt(1.0 - vsqr / clight2);
}


/* ************************************************ */
double Lorentz2Speed(const double lorentz)
/*!
 * Return velocity from Lorentz factor
 *
 ************************************************** */
{
    double clight;
    clight = CONST_c / vn.v_norm;
    return sqrt(1. - 1. / (lorentz * lorentz)) * clight;
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
