/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO main function.
 
  The file main.c contains the PLUTO main function and several other 
  top-level routines.
  main() provides basic code initialization, handles the the principal 
  integration loop and calls the output driver write_data.c.
  Other useful functions contained in this file are Integrate() which does
  the actual integration, GetNextTimeStep() responsible for computing the
  next time step based on the information available at the last time
  level.
  
  We use two slightly different integration loops depending on whether
  asnchrounous I/O has to be performed (macro USE_ASYNC_IO).
  If the macro USE_ASYNC_IO is not defined, the standard 
  integration loop consists of the following steps:
 
   - Check for last step & adjust dt if necessary
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Check output/analysis:  t(n) < tout < t(n)+dt(n)
   - write to disk/call analysis using {U(n), t(n), dt(n)}
   - Advance solution using dt(n): U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Show dominant time step (n)
   - [MPI] Get next time step dt(n+1)
   - [MPI] reduction operations (n)
   - Increment n --> n+1
 
  Otherwise, using Asynchrounous I/O:
 
   - Check for last step & adjust dt
   - check for output/analysis:   t(n) < tout < t(n+1)
     - Write data/call analysis using {U(n), t(n), dt(n)}
   - [MPI] Show dominant time step (n-1)
   - [MPI] reduction operations (n-1)
   - [AIO]: finish writing
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Advance solution using dt(n), U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Get next time step dt(n)
   - Increment n --> n+1

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "globals.h"

#define SHOW_TIME_STEPS  NO   /* -- show time steps due to advection,
                                     diffusion and cooling */

static double GetNextTimeStep (Time_Step *, struct INPUT *, Grid *);
static char *TOTAL_TIME (double dt);
static int Integrate (Data *, Riemann_Solver *, Time_Step *, Grid *);
/* AYW -- 2013-01-08 17:59 JST 
 * Added last_step argument of type char
 * */
static void CheckForOutput (Data *, Input *, Grid *, char);
static void CheckForAnalysis (Data *, Input *, Grid *, char);
/* -- AYW */

/* ********************************************************************* */
int main (int argc, char *argv[])
/*!
 * Start PLUTO, initialize functions, define data structures and 
 * handle the main integration loop.
 *
 * \param [in] argc Argument counts.
 * \param [in] argv Array of pointers to the strings.
 * \return This function return 0 on normal exit.
 *
 *********************************************************************** */
{
  int    nv, idim, err;
  char   first_step=1, last_step=0;
  double scrh;
  Data   data;
  time_t  tbeg, tend;
  /* AYW -- 2012-06-19 10:18 JST
   * Time difference using difftime */
  double t_elapsed;
#if PARALLEL
  double tbeg_mpi, tend_mpi;
#endif
  //struct tm * ptm;
  /* -- AYW */

  Riemann_Solver *Solver;
  Grid      grd[3];
  Time_Step Dts;
  Cmd_Line cmd_line;
  Input  ini;
  Output *output;
/*
  print1 ("sizeof (CMD_LINE)   = %d\n", sizeof(Cmd_Line));
  print1 ("sizeof (DATA)       = %d\n", sizeof(Data));
  print1 ("sizeof (STATE_1D)   = %d\n", sizeof(State_1D));
  print1 ("sizeof (GRID)       = %d\n", sizeof(Grid));
  print1 ("sizeof (TIME_STEP)  = %d\n", sizeof(Time_Step));
  print1 ("sizeof (OUTPUT)     = %d\n", sizeof(Output));
  print1 ("sizeof (INPUT)      = %d\n", sizeof(Input));
  print1 ("sizeof (RUNTIME)    = %d\n", sizeof(Runtime));
  print1 ("sizeof (RGB)        = %d\n", sizeof(RGB));
  print1 ("sizeof (IMAGE)      = %d\n", sizeof(Image));
  print1 ("sizeof (FLOAT_VECT) = %d\n", sizeof(Float_Vect));
  print1 ("sizeof (INDEX)      = %d\n", sizeof(Index));
  print1 ("sizeof (RBOX)       = %d\n", sizeof(RBox));
  QUIT_PLUTO(1);
*/
  #ifdef PARALLEL
   AL_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &prank);
  #endif

  Initialize (argc, argv, &data, &ini, grd, &cmd_line);

/* -- initialize members of Time_Step structure -- */

  Dts.cmax     = ARRAY_1D(NMAX_POINT, double);
  Dts.inv_dta  = 0.0;
  Dts.inv_dtp  = 0.0;
  Dts.dt_cool  = 1.e38;
  Dts.cfl      = ini.cfl;
  Dts.cfl_par  = ini.cfl_par;
  Dts.rmax_par = ini.rmax_par;
  Dts.Nsts     = Dts.Nrkc = 0;
  
  Solver = SetSolver (ini.solv_type);
  

  /* AYW -- 2012-06-26 11:43 JST */
#ifdef PARALLEL  
    if (prank == 0) tbeg_mpi = MPI_Wtime();
    MPI_Bcast(&tbeg_mpi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    time(&tbeg);
#endif
  /* -- AYW */

  /* AYW -- 2012-06-26 11:43 JST */
  //ptm = gmtime ( &tbeg );
  //print1("\n> First Timestamp: %2d:%02d:%02d.\n", (ptm->tm_hour)%24, ptm->tm_min, ptm->tm_sec);
  /* -- AYW */
  g_stepNumber = 0;

/* --------------------------------------------------------
    Check if restart is necessary. 
    If not, write initial condition to disk.
   ------------------------------------------------------- */
   
  if (cmd_line.restart == YES) {
    Restart (&ini, cmd_line.nrestart, DBL_OUTPUT, grd);
  }else if (cmd_line.h5restart == YES){
    Restart (&ini, cmd_line.nrestart, DBL_H5_OUTPUT, grd);
  }else if (cmd_line.write){
    /* AYW -- 2013-01-08 18:05 JST 
     * Arguments modified in CheckFor functions to include last_step*/
    CheckForOutput (&data, &ini, grd, last_step);
    CheckForAnalysis (&data, &ini, grd, last_step);
    //CheckForOutput (&data, &ini, grd);
    //CheckForAnalysis (&data, &ini, grd);
    /* -- AYW */

    #ifdef USE_ASYNC_IO
     Async_EndWriteData (&ini);
    #endif
  }

  print1 ("> Starting computation... \n\n");

/* =====================================================================
          M A I N      L O O P      S T A R T S      H E R E
   ===================================================================== */

#ifndef USE_ASYNC_IO  /* -- Standard loop, don't use Asynchrouns I/O -- */

  while (!last_step){

  /* ------------------------------------------------------
      Check if this is the last integration step:
      - final tstop has been reached: adjust time step 
      - or max number of steps has been reached
     ------------------------------------------------------ */

    if ((g_time + g_dt) >= ini.tstop*(1.0 - 1.e-8)) {
      g_dt   = (ini.tstop - g_time);
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps > 0) {
      last_step = 1;
    }


  /* ------------------------------------------------------
                Dump log information
     ------------------------------------------------------ */

    if (g_stepNumber%ini.log_freq == 0) {
      print1 ("step:%d ; t = %10.4e ; dt = %10.4e ; %d %% ; [%f, %d",
               g_stepNumber, g_time, g_dt, (int)(100.0*g_time/ini.tstop), 
               g_maxMach, g_maxRiemannIter);
      #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
       print1 (", Nsts = %d",Dts.Nsts);
      #endif
      #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
       print1 (", Nrkc = %d",Dts.Nrkc);
      #endif
      print1 ("]\n");      
    }

  /* ------------------------------------------------------
       check if it's time to write or perform analysis
     ------------------------------------------------------ */

    if (!first_step && !last_step && cmd_line.write) {
      /* AYW -- 2013-01-08 18:05 JST 
       * Arguments modified in CheckFor functions to include last_step*/
      CheckForOutput (&data, &ini, grd, last_step);
      CheckForAnalysis (&data, &ini, grd, last_step);
      //CheckForOutput  (&data, &ini, grd);
      //CheckForAnalysis(&data, &ini, grd);
      /* -- AYW */
    }

  /* ------------------------------------------------------
      Advance solution array by a single time step
      g_dt = dt(n)
     ------------------------------------------------------ */

    if (cmd_line.jet != -1) SetJetDomain (&data, cmd_line.jet, ini.log_freq, grd); 
    err = Integrate (&data, Solver, &Dts, grd);
    if (cmd_line.jet != -1) UnsetJetDomain (&data, cmd_line.jet, grd); 

  /* ------------------------------------------------------
       Integration didn't go through. Step must
       be redone from previously saved solution.
     ------------------------------------------------------ */
/*
    if (err != 0){
      print1 ("! Step failed. Re-trying\n");
      zones with problems must be tagged with MINMOD_FLAG and HLL_FLAG
      time step should be halved
      GET_SOL(&data);
    }
*/

  /* ------------------------------------------------------
      Increment time, t(n+1) = t(n) + dt(n)
     ------------------------------------------------------ */

    g_time += g_dt;

  /* ------------------------------------------------------
      Show the time step ratios between the actual g_dt
      and the advection, diffusion and cooling time scales.
     ------------------------------------------------------ */

    #if SHOW_TIME_STEPS == YES
     if (g_stepNumber%ini.log_freq == 0) {
       double cg, dta, dtp, dtc;
       dta = 1.0/Dts.inv_dta;
       dtp = 0.5/Dts.inv_dtp;
       dtc = Dts.dt_cool;
       #ifdef PARALLEL
        MPI_Allreduce (&dta, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dta = cg;

        MPI_Allreduce (&dtp, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtp = cg;

        MPI_Allreduce (&dtc, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtc = cg;
       #endif
       print1 ("[dt/dt(adv) = %10.4e, dt/dt(par) = %10.4e, dt/dt(cool) = %10.4e]\n",
                g_dt/dta, g_dt/dtp, g_dt/dtc);
     }
    #endif


  
  /* ------------------------------------------------------
   * AYW -- 2013-01-08 15:05 JST
   * Check if wallclock time has been reached. Measure 
   * delta t every timestep.
     ------------------------------------------------------ */

#ifdef PARALLEL  
    if (prank == 0) tend_mpi = MPI_Wtime();
    MPI_Bcast(&tend_mpi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    time(&tend);
#endif

#ifdef PARALLEL  
    t_elapsed = tend_mpi - tbeg_mpi;
#else
    t_elapsed = difftime(tend, tbeg);
#endif
    if (t_elapsed >= cmd_line.maxtime && cmd_line.maxtime > 0){
      print1("\n> Max time %f reached.\n", cmd_line.maxtime);
      last_step = 1;
    }
    //else{
    //  ptm = gmtime ( &tend );
    //  print1("\n> Timestamp: %2d:%02d:%02d.\n", (ptm->tm_hour)%24, ptm->tm_min, ptm->tm_sec);
    //  print1("\n> Wall time (s): %f / %f.\n", t_elapsed, cmd_line.maxtime);
    //}
    /* -- AYW */


  /* ------------------------------------------------------
                Get next time step dt(n+1)
     ------------------------------------------------------ */

    g_dt = GetNextTimeStep(&Dts, &ini, grd);

  /* ------------------------------------------------------
          Global MPI reduction operations
     ------------------------------------------------------ */
  
    #ifdef PARALLEL
     MPI_Allreduce (&g_maxMach, &scrh, 1, 
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     g_maxMach = scrh;

     MPI_Allreduce (&g_maxRiemannIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     g_maxRiemannIter = nv;
    #endif

    g_stepNumber++;
    
    first_step = 0;
  }

#else /* Use Asynchrounous I/O */

  while (!last_step){

  /* ------------------------------------------------------
      Check if this is the last integration step:
      - final tstop has been reached: adjust time step 
      - or max number of steps has been reached
     ------------------------------------------------------ */

    if ((g_time + g_dt) >= ini.tstop*(1.0 - 1.e-8)) {
      g_dt   = (ini.tstop - g_time);
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps > 0) {
      last_step = 1;
    }

  /* ------------------------------------------------------
       check if it's time to write or perform analysis
     ------------------------------------------------------ */

    if (!first_step && !last_step && cmd_line.write) {

      /* AYW -- 2013-01-08 18:05 JST 
       * Arguments modified in CheckFor functions to include last_step*/
      CheckForOutput (&data, &ini, grd, last_step);
      CheckForAnalysis (&data, &ini, grd, last_step);
      //CheckForOutput  (&data, &ini, grd);
      //CheckForAnalysis(&data, &ini, grd);
      /* -- AYW */
    }

  /* ------------------------------------------------------
      Show the time step ratios between the actual g_dt
      and the advection, diffusion and cooling time scales.
     ------------------------------------------------------ */

    #if SHOW_TIME_STEPS == YES
     if (!first_step && g_stepNumber%ini.log_freq == 0) {
       double cg, dta, dtp, dtc;
       dta = 1.0/Dts.inv_dta;
       dtp = 0.5/Dts.inv_dtp;
       dtc = Dts.dt_cool;
       #ifdef PARALLEL
        MPI_Allreduce (&dta, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dta = cg;

        MPI_Allreduce (&dtp, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtp = cg;

        MPI_Allreduce (&dtc, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtc = cg;
       #endif
       print1 ("\t[dt/dta = %10.4e, dt/dtp = %10.4e, dt/dtc = %10.4e \n",
                g_dt/dta, g_dt/dtp, g_dt/dtc);
     }
    #endif

  /* ------------------------------------------------------
          Global MPI reduction operations
     ------------------------------------------------------ */
  
    #ifdef PARALLEL
     MPI_Allreduce (&g_maxMach, &scrh, 1, 
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     g_maxMach = scrh;

     MPI_Allreduce (&g_maxRiemannIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     g_maxRiemannIter = nv;
    #endif

  /* ------------------------------------------------------
             Finish writing using Async I/O
     ------------------------------------------------------ */

    #ifdef USE_ASYNC_IO
     Async_EndWriteData (&ini);
    #endif

  /* ------------------------------------------------------
                Dump log information
     ------------------------------------------------------ */

    if (g_stepNumber%ini.log_freq == 0) {
      print1 ("step:%d ; t = %10.4e ; dt = %10.4e ; %d %% ; [%f, %d",
               g_stepNumber, g_time, g_dt, (int)(100.0*g_time/ini.tstop), 
               g_maxMach, g_maxRiemannIter);
      #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
       print1 (", Nsts = %d",Dts.Nsts);
      #endif
      #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
       print1 (", Nrkc = %d",Dts.Nrkc);
      #endif
      print1 ("]\n");      
    }
    
  /* ------------------------------------------------------
      Advance solution array by a single time step
      g_dt = dt(n)
     ------------------------------------------------------ */

    if (cmd_line.jet != -1) SetJetDomain (&data, cmd_line.jet, grd); 
    err = Integrate (&data, Solver, &Dts, grd);
    if (cmd_line.jet != -1) UnsetJetDomain (&data, cmd_line.jet, grd); 

  /* ------------------------------------------------------
       Integration didn't go through. Step must
       be redone from previously saved solution.
     ------------------------------------------------------ */
/*
    if (err != 0){
      print1 ("! Step failed. Re-trying\n");
      zones with problems must be tagged with MINMOD_FLAG and HLL_FLAG
      time step should be halved
      GET_SOL(&data);
    }
*/

  /* ------------------------------------------------------
      Increment time, t(n+1) = t(n) + dt(n)
     ------------------------------------------------------ */

    g_time += g_dt;


  
  /* ------------------------------------------------------
   * AYW -- 2013-01-08 15:05 JST
   * Check if wallclock time has been reached. Measure 
   * delta t every timestep.
     ------------------------------------------------------ */

#ifdef PARALLEL  
    if (prank == 0) tend_mpi = MPI_Wtime();
    MPI_Bcast(&tend_mpi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    time(&tend);
#endif

#ifdef PARALLEL  
    t_elapsed = tend_mpi - tbeg_mpi;
#else
    t_elapsed = difftime(tend, tbeg);
#endif
    if (t_elapsed >= cmd_line.maxtime && cmd_line.maxtime > 0){
      print1("\n> Max time %f reached.\n", cmd_line.maxtime);
      last_step = 1;
    }
    //else{
    //  ptm = gmtime ( &tend );
    //  print1("\n> Timestamp: %2d:%02d:%02d.\n", (ptm->tm_hour)%24, ptm->tm_min, ptm->tm_sec);
    //  print1("\n> Wall time (s): %f / %f.\n", t_elapsed, cmd_line.maxtime);
    //}
    /* -- AYW */


  /* ------------------------------------------------------
                Get next time step dt(n+1)
     ------------------------------------------------------ */

    g_dt = GetNextTimeStep(&Dts, &ini, grd);

    g_stepNumber++;
    
    first_step = 0;
  }
#endif /* USE_ASYNC_IO */

/* =====================================================================
          M A I N       L O O P      E N D S       H E R E 
   ===================================================================== */

  if (cmd_line.write){

    /* AYW -- 2013-01-08 18:05 JST 
     * Arguments modified in CheckFor functions to include last_step*/
    CheckForOutput (&data, &ini, grd, last_step);
    CheckForAnalysis (&data, &ini, grd, last_step);
    //CheckForOutput (&data, &ini, grd);
    //CheckForAnalysis (&data, &ini, grd);
    /* -- AYW */

    #ifdef USE_ASYNC_IO
     Async_EndWriteData (&ini);
    #endif
  }

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   print1  ("\n> Total allocated memory  %6.2f Mb (proc #%d)\n",
             (float)g_usedMem/1.e6,prank);
   MPI_Barrier (MPI_COMM_WORLD);
  #else
   print1  ("\n> Total allocated memory  %6.2f Mb\n",(float)g_usedMem/1.e6);
  #endif

  time(&tend);
  g_dt = difftime(tend, tbeg);
  print1("> Elapsed time             %s\n", TOTAL_TIME (g_dt));
  print1("> Average time/step       %10.2e  (sec)  \n", 
          difftime(tend,tbeg)/(double)g_stepNumber);
  print1("> Local time                %s",asctime(localtime(&tend)));
  print1("> Done\n");

  FreeArray4D ((void *) data.Vc);
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   AL_Finalize ();
  #endif

  return (0);
}
#undef SHOW_TIME_STEPS
/* ******************************************************************** */
int Integrate (Data *d, Riemann_Solver *Solver, Time_Step *Dts, Grid *grid)
/*!
 * Advance equations by a single time-step.

 * \param  d      pointer to PLUTO Data structure;
 * \param  Solver pointer to a Riemann solver function;
 * \param  Dts    pointer to time Step structure;
 * \param  grid   pointer to grid structure.
 * 
 * \return An integer giving success / failure (development).
 * 
 ********************************************************************** */
{
  int idim, err = 0;

  g_maxMach = 0.0;
  g_maxRiemannIter = 0;

/* -------------------------------------------------------
    Initialize max propagation speed in Dedner's approach
   ------------------------------------------------------- */

  #ifdef GLM_MHD  /* -- initialize glm_ch -- */
   GLM_Init (d, Dts, grid);   
   GLM_Source (d->Vc, 0.5*g_dt, grid);
  #endif

  /* ---------------------------------------------
        perform Strang Splitting on direction 
        (if necessary) and sources 
     --------------------------------------------- */

  FlagReset (d);

#ifdef ADD_CYLSOURCE
 CONS_CYLSOLVE(d->Vc, 0.5*g_dt, grid);
#endif

  #ifdef FARGO
   FARGO_ComputeVelocity(d, grid);
  #endif
  if ((g_stepNumber%2) == 0){

    #if DIMENSIONAL_SPLITTING == YES
     for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){
       if (Sweep (d, Solver, Dts, grid) != 0) return (1);
     }
    #else
     if (Unsplit (d, Solver, Dts, grid) != 0) return(1);
    #endif
    SplitSource (d, g_dt, Dts, grid);
  }else{
    SplitSource (d, g_dt, Dts, grid);
    #if DIMENSIONAL_SPLITTING == YES
     for (g_dir = DIMENSIONS - 1; g_dir >= 0; g_dir--){
       if (Sweep (d, Solver, Dts, grid) != 0) return (1);
     }
    #else
     if (Unsplit (d, Solver, Dts, grid) != 0) return(1);
    #endif
  }       


#ifdef ADD_CYLSOURCE
 CONS_CYLSOLVE(d->Vc, 0.5*g_dt, grid);
#endif

  #ifdef GLM_MHD  /* -- GLM source for dt/2 -- */
   GLM_Source (d->Vc, 0.5*g_dt, grid);
  #endif

  return (0); /* -- ok, step achieved -- */
}

/* ******************************************************************** */
char *TOTAL_TIME (double dt)
/*
 * 
 * PURPOSE
 *
 *   convert a floating-point variable (dt, in seconds) to a string 
 *   displaying days:hours:minutes:seconds
 *   
 ********************************************************************** */
{
  static char c[128];
  int days, hours, mins, secs;

  days  = (int) (dt/86400.0);
  hours = (int) ((dt - 86400.0*days)/3600.0);
  mins  = (int) ((dt - 86400.0*days - 3600.0*hours)/60.);
  secs  = (int) (dt - 86400.0*days - 3600.0*hours - 60.0*mins);

  sprintf (c, " %dd:%dh:%dm:%ds", days,hours, mins, secs);
  return (c);
}

/* ********************************************************************* */
double GetNextTimeStep (Time_Step *Dts, struct INPUT *ini, Grid *grid)
/*!
 * Compute and return the time step for the next time level
 * using the information from the previous integration
 * (Dts->inv_dta and Dts->inv_dp).
 *
 * \param [in] Dts pointer to the Time_Step structure
 * \param [in] ini pointer to the Input structure
 * \param [in] grid pointer to array of Grid structures
 *
 * \return The time step for next time level
 *
 *********************************************************************** */
{
  int idim;
  double dt_adv, dt_cool;
  double dtnext, dtnext_glob;
  double dxmin, dtp_dta;
  double dt_par;

  #if    (DIMENSIONS == 1 && INCLUDE_SPLIT_SOURCE == NO) \
      || (DIMENSIONAL_SPLITTING == NO && INCLUDE_SPLIT_SOURCE == NO) 

  #else

  /*  ------------------------------------------------
       For Strang splitting to be 2nd order accurate,  
       change dt only every 2 time steps                 
      ------------------------------------------------ */
 
   if (g_stepNumber%2 == 0) return (g_dt);

  #endif

/* ----------------------------------
        Compute time step
   ---------------------------------- */

  #if (PARABOLIC_FLUX & EXPLICIT)
   dt_adv  = 1.0/(Dts->inv_dta + 2.0*Dts->inv_dtp);
  #else
   dt_adv  = 1.0/Dts->inv_dta;
  #endif
  dt_adv *= ini->cfl;
  dtnext  = dt_adv;

/* --------------------------------------
        Min cell length
   -------------------------------------- */

  dxmin = grid[IDIR].dl_min;
  for (idim = 1; idim < DIMENSIONS; idim++){
    dxmin = MIN(dxmin, grid[idim].dl_min);
  }

/* -------------------------------------------
    Maximum propagation speed for the local
    processor. Global glm_ch will be computed
    later in GLM_Init.
   ------------------------------------------- */

  #ifdef GLM_MHD
   glm_ch = ini->cfl*dxmin/dtnext;
  #endif

/* ------------------------------------------------
    with STS, the ratio between advection (full) 
    and parabolic time steps should not exceed a
    given threshold (200).
   ------------------------------------------------ */
      
  #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
   dt_par  = ini->cfl_par/(2.0*Dts->inv_dtp);
   dtnext *= MIN(1.0, ini->rmax_par/(dt_adv/dt_par));
  #endif

  #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
   dt_par  = ini->cfl_par/(2.0*Dts->inv_dtp);
   dtnext *= MIN(1.0, ini->rmax_par/(dt_adv/dt_par));
  #endif

/* ----------------------------------
      Compute Cooling time step
   ---------------------------------- */

  #if COOLING != NO
   dtnext = MIN(dtnext, Dts->dt_cool);
  #endif

   
/* --------------------------------------------------------------
    allow time step to vary at most by a factor ini->cfl_max_var
   -------------------------------------------------------------- */

  dtnext = MIN(dtnext, ini->cfl_max_var*g_dt);

/* -----------------------------------------------------
      Get min time step among ALL processors
   ----------------------------------------------------- */

  #ifdef PARALLEL
   MPI_Allreduce (&dtnext, &dtnext_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   dtnext = dtnext_glob;
  #endif

/* -----------------------------------------------------
          quit if dt gets too small 
   ----------------------------------------------------- */

  if (dtnext < ini->first_dt*1.e-9){
    print1 ("! dt is too small (%12.6e)!\n", dtnext);
    print1 ("! Cannot continue\n");
    QUIT_PLUTO(1);
  }

/* -----------------------------------------------------
          Reset time step coefficients
   ----------------------------------------------------- */

  DIM_LOOP(idim) Dts->cmax[idim] = 0.0;
  Dts->inv_dta = 0.0;
  Dts->inv_dtp = 0.0;
  Dts->dt_cool = 1.e38;

/* ------------------------------------------------------
     Issue a warning if first_dt has been overestimaed
   ------------------------------------------------------ */

  if (g_stepNumber <= 1 && (ini->first_dt > dtnext/ini->cfl)){
    print1 ("! initial dt exceeds stability limit\n");
  }

  return(dtnext);
}

/* ********************************************************************* */
/* AYW - 2013-01-08 16:56 JST
 * Add last_step as argument
 * */
// void CheckForOutput (Data *d, Input *ini, Grid *grid)
void CheckForOutput (Data *d, Input *ini, Grid *grid, char laststep)
/*!
 *  Check if file output has to be performed.
 *  
 *********************************************************************** */
{
  static int first_call = 1;
  int  n, check_dt, check_dn, check_dclock;
  int  restart_update, last_step;
  double t, tnext;
  Output *output;
  static time_t clock_beg[MAX_OUTPUT_TYPES], clock_end;
  static double tbeg[MAX_OUTPUT_TYPES], tend;
  double dclock;

  restart_update = 0;
  t     = g_time;
  tnext = t + g_dt;
  
  last_step = (fabs(t-ini->tstop) < 1.e-12 ? 1:0);
  /* AYW -- 2013-01-08 18:02 JST 
   * Last step condition from main loop */
  last_step = last_step || laststep;
  /* */

/* -- on first execution initialize
      current beginning time for all output types -- */

  if (first_call){
    #ifdef PARALLEL
     if (prank == 0){
       double tstart;
       tstart = MPI_Wtime();
       for (n = 0; n < MAX_OUTPUT_TYPES; n++) tbeg[n] = tstart;
     }
     MPI_Bcast(tbeg, MAX_OUTPUT_TYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #else
     for (n = 0; n < MAX_OUTPUT_TYPES; n++) time(clock_beg + n);
    #endif
  } 

/* -- get current time -- */

  #ifdef PARALLEL  
   if (prank == 0) tend = MPI_Wtime();
   MPI_Bcast(&tend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #else
   time(&clock_end);
  #endif
  
/* -------------------------------------------------------
          start main loop on outputs
   ------------------------------------------------------- */
   
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    output = ini->output + n;    
    check_dt = check_dn = check_dclock = 0; 

  /* -- check time interval in code units (dt) -- */

    if (output->dt > 0.0){
      check_dt = (int) (tnext/output->dt) - (int)(t/output->dt);
      check_dt = check_dt || g_stepNumber == 0 || last_step;
    }   

  /* -- check time interval in number of steps (dn) -- */

    if (output->dn > 0){
      check_dn = (g_stepNumber%output->dn) == 0;
      check_dn = check_dn || g_stepNumber == 0 || last_step;
    }

  /* -- check time interval in clock time (dclock) -- */

    if (output->dclock > 0.0){
      #ifdef PARALLEL
       dclock = tend - tbeg[n];
      #else
       dclock = difftime(clock_end, clock_beg[n]);
      #endif
      if (dclock >= output->dclock) {
        check_dclock = 1;
        #ifdef PARALLEL
         tbeg[n] = tend;
        #else
         time(clock_beg + n);
        #endif
      }else{ 
        check_dclock = 0;
      }
      check_dclock = check_dclock || g_stepNumber == 0 || last_step;
    }

  /* -- if any of the previous is true dump data to disk -- */

    if (check_dt || check_dn || check_dclock) { 

      #ifdef USE_ASYNC_IO
       if (!strcmp(output->mode,"single_file_async")){
         Async_BegWriteData (d, output, grid);
       }else{
         WriteData(d, output, grid);
       }
      #else     
       WriteData(d, output, grid);
      #endif   

    /* ----------------------------------------------------------
        save the file number of the dbl and dbl.h5 output format
        for writing restart.out once we exit the loop.
       ---------------------------------------------------------- */

      if ((output->type == DBL_OUTPUT) ||
          (output->type == DBL_H5_OUTPUT)) restart_update = 1;
    }
  }

/* -------------------------------------------------------
    Dump restart information if required 

    Note that if both dbl and dbl.h5 formats are used,
    bookkeeping is done using dbl format.
   ------------------------------------------------------- */

  if (restart_update) RestartDump (ini);

  first_call = 0;
}

/* ******************************************************************** */
/* AYW -- 2013-01-08 18:07 JST 
 * Include laststep in argument */
void CheckForAnalysis (Data *d, Input *ini, Grid *grid, char last_step)
//void CheckForAnalysis (Data *d, Input *ini, Grid *grid)
/*
 *
 * PURPOSE 
 *
 *   Check if Analysis needs to be called
 *
 ********************************************************************** */
{
  int check_dt, check_dn;
  double t, tnext;

  t     = g_time;
  tnext = t + g_dt;
  check_dt = (int) (tnext/ini->anl_dt) - (int)(t/ini->anl_dt);
  check_dt = check_dt || g_stepNumber == 0 || fabs(t - ini->tstop) < 1.e-9; 
  check_dt = check_dt && (ini->anl_dt > 0.0);

  check_dn = (g_stepNumber%ini->anl_dn) == 0;
  check_dn = check_dn && (ini->anl_dn > 0);

  /* AYW -- 2013-01-08 18:04 JST 
   * Last step condition from main loop */
  if (check_dt || check_dn || last_step) Analysis (d, grid);
}

