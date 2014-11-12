#include "pluto.h"
#include "pluto_usr.h"
#include "abundances.h"

/* ***************************************************************** */
void PowerLawCooling (Data_Arr VV, real dt, Time_Step *Dts, Grid *grid)
/*
 *
 * PURPOSE:  
 *
 *    Take a source step to account for radiative cooling with
 *    a power-law cooling function.
 *
 *     We directly integrate the following ODE:
 *
 *      dp_cgs/dt_cgs = -(gamma - 1) Lambda(rho_cgs, T_cgs)
 *
 *     where:   Lambda = lambda0 rho_cgs^2/(mu*m_amu,cgs)^2 (T_cgs/T_crit,cgs)^alpha
 *
 *     if the cooling function was a single power law.
 *
 *
 *
 ******************************************************************* */
{
  int   i, j, k;
  real  cost, dE;
  real  rho, p, T, p_f, T_f;
  real  oma, unit_time;

  /* Read cooling parameters file */

  /* Power law cooling mainly valid in regions 10^5
  The values for alpha and lamdbda0 just come from a linear fit to the 
  cooltable.dat data between T = 3*10^4 and 1*10^8 K. See cool_fit_pow.py
  */
  double alpha = -0.37988772;
  double lambda0 = 1.044886865e-22;

  /* Mean molecular mass */
  double vdummy[NVAR];
  real mu = MeanMolecularWeight(vdummy);


  /* A constant that is part of the solution. */
  cost = (g_gamma-1.0)*lambda0/CONST_kB*pow(g_minCoolingTemp, -alpha);

  /* Unit time code to cgs conversion */
  unit_time = g_unitLength/g_unitVelocity;


/*  -------------------------------------------------------------
                Integrate analytically
    -------------------------------------------------------------  */

  //dE = 1.e-18;
  DOM_LOOP(k,j,i){

/*  ----  Find initial temperature in Kelvin  ----  */

    rho = VV[RHO][k][j][i];
    p   = VV[PRS][k][j][i];

    T   = mu*p/rho*KELVIN;

    if (T < g_minCoolingTemp) continue;

/*  ----  Find final temperature  ----  */
     
    /* rho is in code units because the factor in the equation reads 
     * (rho/mu amu), but all other quantities are in cgs */
    oma = 1.-alpha;
    T_f = pow(-cost*oma*rho*dt*unit_time + pow(T, oma), 1./oma);

    if (T_f <= g_minCoolingTemp){ 
      T_f = T - g_maxCoolingRate*(T - g_minCoolingTemp);
    }
    else{
      T_f = T - g_maxCoolingRate*(T - T_f);
    }
   //T_f = MAX(T_f, g_minCoolingTemp);

/*  ----  Update Energy  ----  */

    p_f = T_f*rho/(mu*KELVIN);

    VV[PRS][k][j][i] = p_f;

    //dE = fabs(1.0 - p_f/p) + 1.e-18;
    //Dts->dt_cool = MIN(Dts->dt_cool, dt*g_maxCoolingRate/dE);
  }

}


real MeanMolecularWeight (real *v)
{
  return (0.6156);
}

