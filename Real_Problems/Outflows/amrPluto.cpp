#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "AMR.H"

/* AYW -- 2013-02-06 16:43 JST */
#include "pluto_usr.h"
#include "PLUTOAMR.H"
#include <OldTimer.H>
/* --AYW */

#include "AMRLevelPlutoFactory.H"

#include "UsingNamespace.H"

#ifdef CH_LINUX
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

extern "C" {
#include "globals.h"
}


/* AYW -- 2013-02-20 13:18 JST 
 * Don't know why this is not in the new Chombo anymore */
OldTimer Everything    ("Everything", 0);
OldTimer TimeReadInput ("Read Input",   Everything);
OldTimer TimeSetupAMR  ("Setup AMR",    Everything);
OldTimer TimeRun       ("Run",          Everything);
OldTimer TimeConclude  ("Conclude",     Everything);
/* -- AYW */

#define NOPT      32


/* AYW -- 2013-02-22 16:42 JST */
int maxLevel = 0;
/* -- AYW */

// amrPluto is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrPluto(char *, Cmd_Line *);

/* *********************************************************** */
int main(int a_argc, char* a_argv[])
/*
 *
 *
 *
 *
 *
 *
 *
 ************************************************************* */
{
  Real end_memory, size;
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
//  setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
  prank = rank;
#else
  rank = 0;
  number_procs = 1;
#endif

  if (rank == 0)
  {
    pout() << "Number of procs = " << number_procs << endl;
  }

  /* AYW -- 2013-02-20 13:24 JST
   * Don't know why this isn't in this version of Chombo anymore */
  OldTimer::TimerInit(rank);

  Everything.start();
  /* -- AYW */

  char ini_file[128];
  Cmd_Line cmd_line;

  sprintf (ini_file,"pluto.ini");  // default
  ParseCmdLineArgs(a_argc, a_argv, ini_file, &cmd_line);
  ShowConfig();

#ifdef TRAP_FPE
  enableFpExceptions();
#endif

  // Run amrPluto, i.e., do the computation
  amrPluto(ini_file, &cmd_line);

  /* AYW -- 2013-02-20 13:24 JST
   * Don't know why this isn't in this version of Chombo anymore */
  Everything.stop();
  pout() << endl
         << "Time --- "
         << setw(8) << setprecision(3)
         << setiosflags(ios::fixed)
         << Everything.wc_time()
         << " sec (wall-clock)" << endl << endl;

  /* -- AYW */

#ifdef CH_MPI
  // Exit MPI 
  dumpmemoryatexit();
  MPI_Finalize();
#endif

}

/* *********************************************************** */
void amrPluto(char *ini_file, Cmd_Line *cmd_line)
/*
 *
 *
 *
 *
 *
 *
 ************************************************************* */
{
  Real residentSetSize=0.0;
  Real size=0.0;
  
  Input ini;
  Grid grid[3];
  Setup (&ini, cmd_line, ini_file);
 
  int nghost = GetNghost(NULL);
  for (int idim = 0; idim < DIMENSIONS; idim++) {
    grid[idim].nghost  = nghost;
    grid[idim].np_int  = grid[idim].np_int_glob = ini.npoint[idim];
    grid[idim].np_tot  = grid[idim].np_tot_glob = ini.npoint[idim] + 2*nghost;
    grid[idim].beg     = grid[idim].gbeg = grid[idim].lbeg = nghost;
    grid[idim].end     = grid[idim].gend = grid[idim].lend = (grid[idim].lbeg - 1) + grid[idim].np_int;
    grid[idim].lbound  = ini.lft_bound_side[idim];
    grid[idim].rbound  = ini.rgt_bound_side[idim];
  }
  SetGrid (&ini, grid);
  
  // VERBOSITY
  
  // This determines the amount of diagnositic output generated
  int verbosity = 1; /* -- change through cmd line args ?? -- */
  CH_assert(verbosity >= 0);
  
  // [Grid]
  
  // Set the physical size of the longest dimension of the domain

  Real domainLength = 0.0, xbeg[3], xend[3];

#if (GEOMETRY == SPHERICAL) && (CH_SPACEDIM > 1)
    domainLength = ini.patch_left_node[1][2] - ini.patch_left_node[1][1];
#else
  for (int dir = 0; dir < DIMENSIONS; dir++){
    xbeg[dir] = ini.patch_left_node[dir][1];
    xend[dir] = ini.patch_left_node[dir][2];
    domainLength = MAX(domainLength, xend[dir] - xbeg[dir]); 
  }
#endif
  
  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);

  D_EXPAND(numCells[0] = ini.patch_npoint[IDIR][1]; ,
           numCells[1] = ini.patch_npoint[JDIR][1]; ,
           numCells[2] = ini.patch_npoint[KDIR][1];)

  CH_assert(D_TERM(   (numCells[0] > 0),
                   && (numCells[1] > 0),
                   && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                   && (numCells[1] % 2 == 0),
                   && (numCells[2] % 2 == 0)));

  // REFINEMENT
 
  ParOpen (ini_file);

  //Maximum AMR level limit
  //int maxLevel = 0;
  maxLevel = atoi(ParGet("Levels",1));
  int numReadLevels = MAX(maxLevel, 1);

  // Refinement ratios between levels
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  std::vector<int> refRatios(numReadLevels + 1);
  int totLevels = 1;
  for (int nlev = 0; nlev <= numReadLevels; nlev++){
    refRatios[nlev] = atoi(ParGet("Ref_ratio",nlev+1));
    if (nlev < maxLevel) totLevels *= refRatios[nlev];
  }
 
  pout() << endl;
  pout() << "  Number of levels:      " << maxLevel << endl;
  pout() << "  Equivalent Resolution: " << grid[IDIR].np_int*totLevels;
  D_EXPAND(                                                  , 
           pout() << " x " << grid[JDIR].np_int*totLevels;   ,
           pout() << " x " << grid[KDIR].np_int*totLevels;)
  pout() << endl << endl;

/*
  D_EXPAND(pout() << grid[IDIR].np_int*totLevels; ,
           pout() << " X "+grid[JDIR].np_int*totLevels; ,
           pout() << " X "+grid[KDIR].np_int*totLevels;)
  pout() << endl;
*/
  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals(numReadLevels);
  for (int nlev = 0; nlev < numReadLevels; nlev++)
    regridIntervals[nlev] = atoi(ParGet("Regrid_interval",nlev+1));

  // Threshold that triggers refinement
  Real refineThresh = atof(ParGet("Refine_thresh",1));

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = atoi (ParGet("Tag_buffer_size",1));

  // Minimum dimension of a grid
  int blockFactor = atoi(ParGet("Block_factor",1));

  // Maximum dimension of a grid
  int maxGridSize = atoi(ParGet("Max_grid_size",1));

  NMAX_POINT = maxGridSize + 8;
//  NMAX_POINT = 2*maxGridSize;

  Real fillRatio = atof (ParGet("Fill_ratio",1));

  // [Time] 

  Real cfl         = ini.cfl;         // CFL multiplier
  Real maxDtGrowth = ini.cfl_max_var; // Limit the time step growth
  Real stopTime    = ini.tstop;       // Stop when the simulation time get here
  Real initialDT   = ini.first_dt;    // Initial time step

  Real dtToleranceFactor = 1.1;  // Let the time step grow by this factor 
                                 // above the "maximum" before reducing it
 
  int nstop = 999999;  // Stop after this number of steps
  if (cmd_line->maxsteps > 0) nstop = cmd_line->maxsteps;


  /* AYW -- 2013-02-06 16:39 JST */
  Real maxtime = 1.e30;  // Stop after this number of steps
  if (cmd_line->maxtime > 0) maxtime = cmd_line->maxtime;
  /* -- AYW */

  
  // [Solver]

  Riemann_Solver *Solver;  // SOLVER
  std::string riem = ini.solv_type;
  Solver = SetSolver(riem.c_str());
 
  Real fixedDt = -1; // Determine if a fixed or variable time step will be used

  // [Boundary] 

  int  leftbound[3];
  int  rightbound[3];
  vector<int> isPeriodica(SpaceDim,0);

  leftbound[0]  = ini.lft_bound_side[0];
  leftbound[1]  = ini.lft_bound_side[1];
  leftbound[2]  = ini.lft_bound_side[2];
  rightbound[0] = ini.rgt_bound_side[0];
  rightbound[1] = ini.rgt_bound_side[1];
  rightbound[2] = ini.rgt_bound_side[2];

  for (int dir = 0; dir < DIMENSIONS; dir++){
    if (leftbound[dir] == PERIODIC || rightbound[dir] == PERIODIC){
      isPeriodica[dir] = 1;
    } else {
      isPeriodica[dir] = 0;
    }
  }

  bool isPeriodic[SpaceDim];
  // convert periodic from int->bool
  for (int dim = 0; dim < SpaceDim; dim++)
     {
       isPeriodic[dim] = (isPeriodica[dim] == 1);
       if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
         pout() << "Using Periodic BCs in direction: " << dim << endl;
     }

  // [Output]

  // Set up checkpointing
  Real checkpointPeriod = atof(ParGet("Checkpoint_interval",1));
  int checkpointInterval = atoi(ParGet("Checkpoint_interval",2));
 
  // Set up plot file writing
  Real plotPeriod = atof(ParGet("Plot_interval",1));
  int plotInterval = atoi(ParGet("Plot_interval",2));

  if (cmd_line->write == NO) {
    plotPeriod = checkpointPeriod = -1.0;
    plotInterval = checkpointInterval = -1;
  }
 
  // Set up output files ---> after definition of amr
 
  // Restart ---> after definition of amr
 

  // [Parameters]

  for (int i = 0; i < USER_DEF_PARAMETERS; ++i) g_inputParam[i] = ini.aux[i];

  // Call Init once so all processors share global variables
  // assigniments such as g_gamma, g_unitDensity, and so on.
  // This is necessary since not all processors call Init
  // from Startup.

  double u[256];
  Init (u, grid[IDIR].x[0], grid[JDIR].x[0], grid[KDIR].x[0]);

  // Print the parameters
  if ( verbosity >= 2 )
    {
      pout() << "maximum step = " << nstop << endl;
      pout() << "maximum time = " << stopTime << endl;

      pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                               numCells[1] << "  " <<,
                                               numCells[2] << ) endl;

      pout() << "maximum level = " << maxLevel << endl;

      pout() << "refinement ratio = ";
      for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
      pout() << endl;

      pout() << "regrid interval = ";
      for (int i = 0; i < regridIntervals.size(); ++i) pout() << regridIntervals[i] << " ";
      pout() << endl;

      pout() << "refinement threshold = " << refineThresh << endl;

      pout() << "blocking factor = " << blockFactor << endl;
      pout() << "max grid size = " << maxGridSize << endl;
      pout() << "fill ratio = " << fillRatio << endl;

      pout() << "checkpoint interval = " << checkpointInterval << endl;
      pout() << "plot interval = " << plotInterval << endl;
      pout() << "CFL = " << cfl << endl;
      pout() << "initial dt = " << initialDT << endl;
      if (fixedDt > 0)
        {
          pout() << "fixed dt = " << fixedDt << endl;
        }
      pout() << "maximum dt growth = " << maxDtGrowth << endl;
      pout() << "dt tolerance factor = " << dtToleranceFactor << endl;
    }

  ShowUnits();

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  initialDT = initialDT*probDomain.domainBox().longside();

  PatchPluto* patchPluto =  static_cast<PatchPluto*>(new PatchPluto());

  // Set up the patch integrator
  patchPluto->setBoundary(leftbound,
                          rightbound);
  patchPluto->setRiemann(Solver);

  // Set up the AMRLevel... factory
  AMRLevelPlutoFactory amrGodFact;

  amrGodFact.define(cfl,
                    domainLength,
                    verbosity,
                    refineThresh,
                    tagBufferSize,
                    initialDT,
                    patchPluto);

  /* AYW -- 2013-02-06 16:43 JST */
  //AMR amr;
  PLUTOAMR amr;
  /* -- AYW */

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrGodFact);

  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // Set output parameters
  amr.checkpointPeriod(checkpointPeriod);
  amr.checkpointInterval(checkpointInterval);
  amr.plotPeriod(plotPeriod);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);

  // Set up output files
  std::string prefix = "data.";
  amr.plotPrefix(prefix);

  // Set up checkpoint files
  prefix = "chk.";
  amr.checkpointPrefix(prefix);
  
  amr.verbosity(verbosity);

  // Set up input files

  if (cmd_line->restart == NO && cmd_line->h5restart == NO){
    if (!ParQuery("fixed_hierarchy")) {

      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    } else       {
      MayDay::Error("Fixed_hierarchy option still disabled");
    }
  } else {
    char restartFile[128];
    sprintf (restartFile,"chk.%04d.hdf5",cmd_line->nrestart);

    pout() << endl << "> Restarting from file " << restartFile << endl;

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
    // Create dummies to be passed to startup
    // startup should be called at restart too
/*
    Box dummyBox (IntVect::Zero,
                  IntVect(D_DECL(NMAX_POINT-1,
                                 NMAX_POINT-1,
                                 NMAX_POINT-1)));
    FArrayBox dummyData(dummyBox,NVAR);
    dummyData.setVal(0.0);
    patchPluto->startup(dummyData);
*/
#else
    MayDay::Error("amrPluto restart only defined with hdf5");
#endif
  }

  // Run the computation
  pout() << "> Starting Computation" << endl;
  /* AYW -- 2013-02-06 16:42 JST */
  Real cur_walltime = Everything.wc_time();
  amr.run(stopTime, nstop, cur_walltime, maxtime);
  //amr.run(stopTime,nstop);
  /* -- AYW */

  pout() << "> Done\n";

  // Output the last plot file and statistics
  amr.conclude();

}


/* ******************************************************************* */
void print (const char *fmt, ...)
/* 
 * 
 * General purpose print function used by PLUTO-Chombo.
 * The corresponding static grid version is given in tools.c
 *
 ********************************************************************* */
{
  char buffer[256];
  va_list args;
  va_start (args, fmt);
  vsprintf (buffer,fmt, args);
  pout() << buffer;
}
void print1 (const char *fmt, ...)
{
  char buffer[256];

  if (prank != 0) return;
  va_list args;
  va_start (args, fmt);
  vsprintf (buffer,fmt, args);
  pout() << buffer;
}

#undef NOPT
