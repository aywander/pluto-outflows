#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/* AYW -- 2013-02-06 16:43 JST
 * This is to ensure that this routine is only compiled if Chombo is being used.
 * I chose a pretty random Chombo PP variable for this check.*/
#ifdef CH_USE_64
/* -- AYW */


/* AYW -- 2013-02-20 12:34 JST */
//#ifndef _AMR_H_
//#define _AMR_H_
#ifndef _PLUTOAMR_H_
#define _PLUTOAMR_H_
/* -- AYW */

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
#include "Scheduler.H"

/* AYW -- 2013-02-06 16:43 JST */
#include "AMR.H"
/* -- AYW */

#include "NamespaceHeader.H"


/// Framework for Berger-Oliger timestepping for AMR
/**
   This class is a framework for Berger-Oliger timestepping for
   adaptive mesh refinement of time-dependent problems.  It is
   applicable to both hyperbolic and parabolic problems.  It
   represents a hierarchy of levels of refinement as a collection of
   AMRLevel objects.

   The usage pattern is this:

   <ol>
     <li>Call define() to define the stuff that does not change throughout
         the run (maxlevel, refinement ratios, domain, and operator).
     <li>Modify any parameters you like (blocking factor, max grid size, etc.)
         using parameter-setting functions.
     <li>Call any one of the three setup functions (setupForRestart(),
         setupForNewAMRRun(), or setupForFixedHierarchyRun())
         so AMR can set up all its internal data structures.
     <li>Call run() to run the calculation.
     <li>Call conclude() to output (to stdout) how many cells were updated
         and so forth.
   </ol>
*/
class PLUTOAMR : public AMR
{
public:

  using AMR::run;
  void run(Real a_max_time, int a_max_step) { AMR::run(a_max_time, a_max_step); }
  void run(Real a_max_time, int a_max_step, Real a_start_wtime, float a_max_wtime);

};

#include "NamespaceFooter.H"
#endif

/* AYW -- 2013-02-06 16:43 JST */
#endif
/* -- AYW */
