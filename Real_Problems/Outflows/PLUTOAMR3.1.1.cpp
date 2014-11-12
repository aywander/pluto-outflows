#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/* AYW -- 2012-06-25 15:12 JST 
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

#include "REAL.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "BRMeshRefine.H"
#include "DisjointBoxLayout.H"
#include "CH_HDF5.H"
#include "AMRLevel.H"
#include "parstream.H"
#include "Tuple.H"
#include "BoxIterator.H"
//#include "AMR.H"
#include "CH_Timer.H"

/* AYW -- 2012-06-19 13:20 JST */
#include "PLUTOAMR.H"
#include "pluto_usr.h"
/* --AYW */

#include "NamespaceHeader.H"

#ifdef CH_USE_TIMER
using namespace Chombo;
#endif

#define SMALL_TIME 1.0e-8


/* AYW -- 2012-06-19 13:45 JST */
extern OldTimer Everything;
/* --AYW */


void PLUTOAMR::run(Real a_max_time, int a_max_step, Real a_start_wtime, float a_max_wtime)
{
  CH_TIME("AMR::run");

  CH_assert(isDefined());
  CH_assert(isSetUp());

  /* AYW - just testing: Remove this later */ 
#ifdef CH_USE_TIMER
  double last_timestep_time = 0 ;
#endif

  /* AYW -- 2012-06-19 14:56 JST */
  Everything.stop();
  double cur_walltime = Everything.wc_time();
  Everything.start();
  /* -- AYW */

  if (m_verbosity >= 3)
  {
    pout() << "AMR::coarseTimeStep:" << endl;
    pout() << "max_time = " << a_max_time << endl;
    pout() << "max_step = " << a_max_step << endl;
    /* AYW -- 2012-06-19 11:47 JST */
    pout() << "start wall time = " << a_start_wtime << endl;
    pout() << "max wall time = " << a_max_wtime << endl;
    pout() << "current wall time = " << cur_walltime << endl;
    /* --AYW */
  }

  Real old_dt_base = m_dt_base;

  for ( ; (m_cur_step < a_max_step) &&
      /* AYW -- 2012-06-19 11:49 JST */
      (cur_walltime < a_max_wtime) &&
      //(a_start_wtime + last_timestep_time < a_max_wtime) &&
      /* --AYW */
      (a_max_time - m_cur_time > m_time_eps*m_dt_base);
      ++m_cur_step, m_cur_time += old_dt_base)
  {
#ifdef CH_USE_TIMER
    m_timer->start() ;
#endif
    old_dt_base = m_dt_base;
    for (int level = 0; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->time(m_cur_time);
    }

    if (m_checkpoint_time > 0.0) {
      Real a_next_time = m_cur_time + old_dt_base;
      int check_dt;

      if (m_cur_step > 0) check_dt = (int) (a_next_time/m_checkpoint_time) - (int)(m_cur_time/m_checkpoint_time);
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

    if ((m_checkpoint_interval > 0)      &&
        (m_lastcheck_step != m_cur_step) &&
        (m_restart_step != m_cur_step)   &&
        (m_cur_step % m_checkpoint_interval == 0))
    {
      m_num_check += 1;
      writeCheckpointFile();
      m_lastcheck_step= m_cur_step;
    }

    if (m_plot_time > 0.0) {
      Real a_next_time = m_cur_time + old_dt_base;
      int check_dt;

      if (m_cur_step > 0) check_dt = (int) (a_next_time/m_plot_time) - (int)(m_cur_time/m_plot_time);
      else check_dt = 1;

      if (check_dt == 1) {
        m_num_plot += 1;
        writePlotFile();
      }
    }

    if ((m_plot_interval > 0) &&
        (m_cur_step % m_plot_interval == 0))
    {
      m_num_plot += 1;
      writePlotFile();
    }

    int level = 0;
    int stepsLeft = 0;
    bool timeBoundary = true;
    (void)timeStep(level,stepsLeft,timeBoundary);

    assignDt();
#ifdef CH_USE_TIMER
    m_timer->stop();
#endif
    if(m_verbosity >= 1)
    {
      std::ios::fmtflags origFlags = pout().flags();
      int origWidth = pout().width();
      int origPrecision = pout().precision();
      pout() << resetiosflags(ios::fixed)
        << "coarse time step " << setw(3) << m_cur_step
        << "  old time = " << setw(12)
        << setprecision(6)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << m_cur_time
        << "  old dt = "   << setw(12)
        << setprecision(6)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << old_dt_base
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
#ifdef CH_USE_TIMER
    last_timestep_time = m_timer->wc_time() ;
#endif
    /* AYW -- 2012-06-19 14:56 JST */
    Everything.stop();
    cur_walltime = Everything.wc_time();
    Everything.start();

    /* AYW -- 2012-06-19 11:47 JST */
    //pout() << "start wall time = " << a_start_wtime << endl;
    //pout() << "max wall time = " << a_max_wtime << endl;
    //pout() << "current wall time = " << cur_walltime << endl;
    /* --AYW */

    /* -- AYW */
  }
}

#ifdef CH_USE_TIMER
Timer * PLUTOAMR::timer(Timer *a_timer )
{
  Timer * old_timer = m_timer ;
  if( a_timer != NULL ){
    m_timer = a_timer ;
  }
  return old_timer ;
}
#endif

#include "NamespaceFooter.H"

/* AYW -- 2012-06-25 15:12 JST */
#endif
/* -- AYW */
