#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/* AYW -- 2013-02-06 17:18 JST
 * This is to ensure that this routine is only compiled if Chombo is being used.
 * I chose a pretty random Chombo PP variable for this check.*/
#ifdef CH_USE_64
/* -- AYW */


#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>

#ifndef CH_DISABLE_SIGNALS
#include <signal.h>  // For handling Ctrl-C.
#endif

#include "IntVectSet.H"
#include "parstream.H"
#include "Tuple.H"
#include "BoxIterator.H"
#include "AMR.H"
#include "CH_Timer.H"

/* AYW -- 2013-02-20 13:41 JST */
#include "PLUTOAMR.H"
#include <OldTimer.H>
/* -- AYW */

#include "NamespaceHeader.H"

#ifdef CH_USE_TIMER
using namespace Chombo;
#endif

#define SMALL_TIME 1.0e-8


//-----------------------------------------------------------------------
// RUDIMENTARY SIGNAL HANDLING -- when control-C is pressed, what happens?
// If we're taking timesteps, we don't want to drop everything as soon as
// we hear a Ctrl-C. This signal handler allows us to set a flag to determine
// whether someone has pressed Ctrl-C and then act on it when we please.
// For the record, I'm putting this in here in order to be able to break out
// of a run within a Python interpreter. Python sets up its own SIGINT handler,
// so this enables us to get it back and treat it safely in the context of
// AMR calculations. -Jeff Johnson, 5/19/2010.

// This flag is set to true if an interrupt signal was intercepted, false otherwise.
static bool s_interrupted = false;

#ifndef CH_DISABLE_SIGNALS
//-----------------------------------------------------------------------
static void
handleCtrlC(int signum)
{
  s_interrupted = true;
}
//-----------------------------------------------------------------------
#endif

extern Real ch_max, ch_max_loc, coeff_dl_min;
extern int  use_glm;

/* AYW -- 2013-02-06 16:44 JST */
extern OldTimer Everything;
/* --AYW */

void PLUTOAMR::run(Real a_max_time, int a_max_step, Real a_start_wtime, float a_max_wtime)
{
  CH_TIME("PLUTOAMR::run");

  CH_assert(isDefined());
  CH_assert(isSetUp());

#ifdef CH_USE_TIMER
  double last_timestep_time = 0 ;
#endif

  /* AYW -- 2013-02-06 16:44 JST */
  //Everything.stop();
  double cur_walltime = Everything.wc_time();
  //Everything.start();
  /* -- AYW */

  if (m_verbosity >= 3)
    {
      pout() << "AMR::coarseTimeStep:" << endl;
      pout() << "max_time = " << a_max_time << endl;
      pout() << "max_step = " << a_max_step << endl;
    /* AYW -- 2013-02-06 16:44 JST */
    pout() << "start wall time = " << a_start_wtime << endl;
    pout() << "max wall time = " << a_max_wtime << endl;
    pout() << "current wall time = " << cur_walltime << endl;
    /* --AYW */

    }

#ifndef CH_DISABLE_SIGNALS
  // Replace the current signal handler for Ctrl-C for the duration of
  // this method.
  struct sigaction ctrlC, oldCtrlC;
  ctrlC.sa_handler = handleCtrlC;
  sigemptyset(&ctrlC.sa_mask);
  ctrlC.sa_flags = 0;
  sigaction(SIGINT, &ctrlC, &oldCtrlC);
#endif

  Real old_dt_base = m_dt_base;
  int  last_timestep = 0; /* 0  = not last time step
                             1  = last time step
                            -1  = quit           */

  if (m_cur_step == 0) ch_max = 1.e-10;

  for ( ; (m_cur_step < a_max_step) && 
      /* AYW -- 2013-02-06 16:44 JST */
      (cur_walltime < a_max_wtime) &&
      //(a_start_wtime + last_timestep_time < a_max_wtime) &&
      /* --AYW */
      (last_timestep != -1);
        ++m_cur_step, m_cur_time += old_dt_base)
      {
        s_step = m_cur_step;
#ifdef CH_USE_TIMER
      m_timer->start() ;
#endif
      old_dt_base = m_dt_base;
      for (int level = 0; level <= m_max_level; ++level)
        {
          m_amrlevels[level]->time(m_cur_time);
        }

       Real a_next_time = m_cur_time + old_dt_base;
      // Drop a checkpoint file if we've gone enough steps. 
      if ((m_checkpoint_interval > 0)      &&
          (m_lastcheck_step != m_cur_step) &&
          (m_restart_step != m_cur_step)   &&
          (m_cur_step % m_checkpoint_interval == 0))
        {
          m_num_check += 1;
          writeCheckpointFile();
          m_lastcheck_step= m_cur_step;
        }

      // Drop a checkpoint file if enough time has passed. 
      if (m_checkpoint_period > 0.0) {
       int check_dt;

       if (m_cur_step > 0) check_dt = (int) (a_next_time/m_checkpoint_period) - (int)(m_cur_time/m_checkpoint_period);
       else check_dt = 1;

      if ((m_lastcheck_step != m_cur_step) &&
          (m_restart_step != m_cur_step)   &&
          (check_dt == 1))
        {
          m_num_check += 1;
          writeCheckpointFile();
          m_lastcheck_step= m_cur_step;
        }
       }

      // Plot if we've gone enough steps.
      if ((m_plot_interval > 0) &&
          (m_cur_step % m_plot_interval == 0))
        {
          m_num_plot += 1;
          writePlotFile();
        }

      // Plot if enough time has passed.
      if (m_plot_period > 0.0) {
       int check_dt;

       if (m_cur_step > 0) check_dt = (int) (a_next_time/m_plot_period) - (int)(m_cur_time/m_plot_period);
       else check_dt = 1;

       if (check_dt == 1) {
         m_num_plot += 1;
         writePlotFile();
       }
      }

      if (m_verbosity >= 1)
        {
          std::ios::fmtflags origFlags = pout().flags();
          int origWidth = pout().width();
          int origPrecision = pout().precision();
          pout() << resetiosflags(ios::fixed)
                 << "coarse step: " << setw(3) << m_cur_step
                 << ";  t = " << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << m_cur_time
                 << ";  dt = "   << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << old_dt_base << ";"
#ifdef CH_USE_TIMER
                 << "  wallclocktime = "
                 << m_timer->wc_time() - last_timestep_time
#endif
                 << resetiosflags(ios::scientific)
                 << endl;
          pout().flags(origFlags);
          pout().width(origWidth);
          pout().precision(origPrecision);
        }

      // Call any scheduled functions. This is placed here so that
      // the plotter function can assume plot files have already been dumped.
      if (!m_scheduler.isNull())
        m_scheduler->execute(m_cur_step, m_cur_time);

      ch_max_loc = 0.0;

      int level = 0;
      int stepsLeft = 0;
      bool timeBoundary = true;
      (void)timeStep(level,stepsLeft,timeBoundary);

      if (use_glm == 1) {
       #ifdef CH_MPI
        int result = MPI_Allreduce(&ch_max_loc, &ch_max, 1, MPI_CH_REAL,
                                   MPI_MAX, Chombo_MPI::comm);
        if(result != MPI_SUCCESS){ //bark!!!
          MayDay::Error("sorry, but I had a communcation error on Ch");
        }
       #else
        ch_max = ch_max_loc;
       #endif
      }

      if (last_timestep == 1) last_timestep = -1;
      else{
        if ((a_next_time + m_dt_new[0]) > a_max_time){
          m_dt_new[0] = a_max_time - a_next_time;
          last_timestep = 1;
        }
      }

      assignDt();
#ifdef CH_USE_TIMER
      m_timer->stop();
#endif

#ifdef CH_USE_TIMER
      last_timestep_time = m_timer->wc_time() ;
#endif

      // If we have assigned a signal handler for interrupts, check for
      // an interrupt and call the handler.
      if (s_interrupted)
        break;

    /* AYW -- 2013-02-06 16:44 JST */
    //Everything.stop();
    cur_walltime = Everything.wc_time();
    //Everything.start();

    /* AYW -- 2013-02-06 16:44 JST */
    //pout() << "start wall time = " << a_start_wtime << endl;
    //pout() << "max wall time = " << a_max_wtime << endl;
    //pout() << "current wall time = " << cur_walltime << endl;
    /* --AYW */

    /* -- AYW */

    }

#ifndef CH_DISABLE_SIGNALS
  // Re-instate the old Ctrl-C handler.
  sigaction(SIGINT, &oldCtrlC, NULL);

  // If Ctrl-C was uttered, clear the interrupt flag and notify the old
  // handler.
  if (s_interrupted)
  {
    s_interrupted = false;
    oldCtrlC.sa_handler(SIGINT);
  }
#endif
}


#include "NamespaceFooter.H"

/* AYW -- 2013-02-06 16:44 JST */
#endif
/* -- AYW */

