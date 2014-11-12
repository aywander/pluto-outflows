#include "pluto.h"
#include "pluto_usr.h"
#include "abundances.h"
#include "interpolation.h"

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
  real  rho, p, T, p_f, T_f;
  //real  dE;

  /* Read cooling parameters file */
  FILE *fcool;
  static real ivls[COOLING_NFIT+1], fpar[COOLING_NFIT][COOLING_OFIT+1];
  static int read = 0;
  int it, ot;
  if (read == 0){
    print1 (" > Reading cooling function fit parameters from file...\n");
    fcool = fopen(COOL_FIT_FILE, "r"); 
    if (fcool == NULL){
      print1 ("! %s does not exists\n", COOL_FIT_FILE);
      QUIT_PLUTO(1);
    }
    for (it=0; it<COOLING_NFIT+1; ++it){
      fscanf(fcool, "%lf", &ivls[it]);
    }
    for (it=0; it<COOLING_NFIT; ++it){
      for (ot=0; ot<COOLING_OFIT+1; ++ot){
        fscanf(fcool, "%lf", &fpar[it][ot]);
      }
    }
    read = 1;

    /* Replace lowest temperature in interval with the lower cooling limit */
    ivls[0] = g_minCoolingTemp;
  }

  /* Mean molecular mass */
  double vdummy[NVAR];
  real mu = MeanMolecularWeight(vdummy);

  /* A constant that is part of the solution. */

  /* Unit time code to cgs conversion */
  real unit_time = g_unitLength/g_unitVelocity;

  /* Calculate density, temperature, and time independent 
   * paramters for all fit intervals here first.
   * Find integration time steps across intervals. 
   * This is density dependent, but time independent. 
   * It only depends on the temperatures either side 
   * of the intervals. Hence, we calculate the deltat 
   * per unit density. deltat is in cgs units. */
  double alpha[COOLING_NFIT], lambda0[COOLING_NFIT];
  double deltat[COOLING_NFIT];
  double oma[COOLING_NFIT], cost[COOLING_NFIT];
  for (it=0; it<COOLING_NFIT; ++it){
    alpha[it] = fpar[it][0];
    lambda0[it] = pow(10, fpar[it][1]);
    oma[it] = 1.-alpha[it];
    cost[it] = (g_gamma-1.0)*lambda0[it]/CONST_kB*pow(g_minCoolingTemp, -alpha[it]);
    deltat[it] = (pow(ivls[it], oma[it]) - pow(ivls[it+1], oma[it]))/(-cost[it]*oma[it]);
  }



  /*  -------------------------------------------------------------
      Integrate analytically
      -------------------------------------------------------------  */


  /* A time interval summation variable */
  double sdt;
  int it0;

  DOM_LOOP(k,j,i){


    /*  ----  Find initial temperature in Kelvin  ----  */

    rho = VV[RHO][k][j][i];
    p   = VV[PRS][k][j][i];

    T   = mu*p/rho*KELVIN;

    if (T > g_minCoolingTemp){

      /* Find interval in which current T lies, and calculate
       * integration time interval to the lower T boundary of 
       * that interval.
       *    rho is in code units because the factor in the equation
       * reads (rho/mu amu), but all other quantities are in cgs */
      it0 = it = hunter(ivls, COOLING_NFIT+1, T);
      sdt = (pow(ivls[it], oma[it]) - pow(T, oma[it]))/(-cost[it]*oma[it]*rho);

      /* Loop over fitting intervals. Keep going until we cross
       * full integration time interval. */
      while (sdt < dt*unit_time && it > 0){
        it -= 1;
        sdt += deltat[it]/rho;
      }

      /*  ----  Find final temperature  ----  */
      /* rho is in code units because the factor in the equation reads 
       * (rho/mu amu), but all other quantities are in cgs.
       * Also, go back one interval first in the case it!=it0. */
      if (it == it0){ 
        T_f = pow(MAX(-cost[it]*oma[it]*rho*dt*unit_time + 
              pow(T, oma[it]), 0), 1./oma[it]);
      }
      else {
        sdt -= deltat[it]/rho;
        T_f = pow(MAX(-cost[it]*oma[it]*rho*(dt*unit_time - sdt) + 
              pow(ivls[it+1], oma[it]), 0), 1./oma[it]);
      }


      /* Limit cooling -- this is a bit of a fudge */
      if (T_f <= g_minCoolingTemp){ 
        T_f = T - g_maxCoolingRate*(T - g_minCoolingTemp);
      }
      else{
        T_f = T - g_maxCoolingRate*(T - T_f);
      }


      /*  ----  Update Energy  ----  */

      p_f = T_f*rho/(mu*KELVIN);

      VV[PRS][k][j][i] = p_f;

      //Dts->dt_cool = MIN(Dts->dt_cool, dt*g_maxCoolingRate/dE);
    }
  }

}


real MeanMolecularWeight (real *v)
{
  return (0.6156);
}

