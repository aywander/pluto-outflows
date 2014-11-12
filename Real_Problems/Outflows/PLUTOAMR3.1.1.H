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


#ifndef _PLUTOAMR_H_
#define _PLUTOAMR_H_

#include <iostream>
#include <string>
#include <ctime>

#include "REAL.H"
#include "Vector.H"
#include "AMRLevel.H"
#include "AMRLevelFactory.H"
#include "BRMeshRefine.H"
#include "ProblemDomain.H"
#include "Box.H"
#include "CH_HDF5.H"
#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

/* AYW -- 2012-06-19 14:14 JST */
#include "AMR.H"
/* -- AYW */

#include "NamespaceHeader.H"

/// Framework for Berger-Oliger timestepping for AMR
//
// AYW -- 
// This is is a slightly modified version for use with PLUTO. 
// The AMR class doesn't implement virtual functions but we override the run function
// anyway. 
//
// Things changed are:
//  - call to AMR.run whose argument now includes
//    a_offset_wtime and a_max_wtime which are the offset time AMR::run was called
//    job (so that we can take into account initalization time etc) and the maximum 
//    wall clock time.
//  - the integration loop in run checks current time versus finish time
//
//  -- AYW
/**
   This class is a framework for Berger-Oliger timestepping for
   adaptive mesh refinement of time-dependent problems.  It is
   applicable to both hyperbolic and parabolic problems.  It
   represents a hierarchy of levels of refinement as a collection of
   AMRLevel objects.

   The usage pattern is this:

   <ol>
     <li>Call define to define the stuff that does not change throughout
         the run (maxlevel, refinement ratios, domain, and operator).
     <li>Modify any parameters you like (blocking factor and so forth
         using access functions.
     <li>Call any one of the three setup functions so AMR can set up all its
         internal data structures.
     <li>Call run to run the calculation.
     <li>Call conclude to output (to stdout) how many cells were updated
         and so forth.
   </ol>
*/
class PLUTOAMR : public AMR
{
public:


  /**
     Runs the calculation.  User must have previously called
     both the define function and a setup function.
     The function AMR::run will be hidden so need to unhide it 
     with using AMR::run first in order to overload.
  */
  using AMR::run;
  void run(Real a_max_time, int a_max_step) { AMR::run(a_max_time, a_max_step); }
  void run(Real a_max_time, int a_max_step, Real a_start_wtime, float a_max_wtime);

#ifdef CH_USE_TIMER
  Chombo::Timer *m_timer;  //assumes the application manages the memory
#endif

};

#include "NamespaceFooter.H"
#endif

/* AYW -- 2012-06-25 15:12 JST */
#endif
/* -- AYW */
