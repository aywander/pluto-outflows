#!/bin/bash
#PBS -P n72
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l mem=960GB
#PBS -l ncpus=240
#PBS -lstorage=scratch/n72+gdata/n72
#PBS -r y
#PBS -l wd
#PBS -m abe
#PBS -M alexander.y.wagner@gmail.com

#_______________________________________________________
# Gadi limits (48 cores with 192 GB memory per node)
#
# 48 hours for 1-255 cores
# 24 hours for 256-511 cores
# 10 hours for 512-1023 cores
# 5 hours for 1024-56960 cores
#_______________________________________________________
# Some user settings

# TODO: Refactor some of the tasks to functions. This is to simplify adoption of this submission script for other programs.

# PLUTO arguments
PLUTOARGS=()

# The value to the pluto -maxsteps option.
# Leave empty if this option is not desired.
MAXSTEPS=

# The value to the pluto -maxtime option, which is the walltime (sec)
# at which PLUTO quits its integration loop (not available in the
# original PLUTO release). Leave empty if this option is
# not desired. Currently, a simulation ending due to the MAXTIME condition
# does not guarantee output of a restart file (to be amended)
# It may be more desirable to set the wall clock output rate of of dbls
# to ensure that restart files are dumped (which is available since PLUTO 4.3)
MAXTIME="1.55h"

# The job number to start at. If START_NJOB > 1
# then submission is considered a "manual restart"
# and the script goes into the restart logic, and
# starts pluto with -restart. Default is 1.
START_NJOB=1

# File number to restart from. Leave empty if the
# last file number is to be used.
# Specify the initial as %4d, e.g. 0023, or just 23.
RESTART_FNUM=

# Target job number. Script auto resubmits until
# the job counter (NJOB) reaches this number.
# Note, NJOBS is *not* necessarily the number of jobs
# to be (re)submitted, but it is the number up to which
# the job counter NJOB should go, starting at START_NJOB.
# Therefore, we require NJOBS > START_NJOB.
NJOBS=10

# Maximum number of crashes allowed.
NCRASHMAX=3

# PLUTO ini file used
inifile="pluto.ini"

# Working directory (reference directory for this script)
WDIR="."

# PLUTO Executable. If left blank, defaults to "./pluto".
EXEC=

# Load modules
module load openmpi
module load hdf5/1.10.5p

# Which mpirun to use. If left blank uses mpirun in path.
MPIRUN=
#module load hpctoolkit
#export HPCRUN_EVENT_LIST="WALLCLOCK@500"
#MPIRUN="mpirun -np 128 hpcrun"

# Grid mode: Chombo or empty for native PLUTO grid
GRID="native"


#_________________________________________________
# Settings that probably don't need changing

# Some Chombo or native grid dependent settings.
#
#   Checkpointfile strings/patterns. Checkpoint files are
# searched by $chksfx*$chkpfx. The pattern $chknumpat is
# extracted from the file name to determine restart #.
# Note that output in dbl format can be single or multiple files
# Change the chkpfx accordingly: "data." for single file, "rho."
# for multiple files.
#
#   goodmsg: Message in output file that confirms run ended fine
# The first below is for native grid runs, the second for
# Chombo runs. If no such check is desired, leave empty.
#
#   logfile: The log file to which output is written to
#
if [[ "$GRID" == "native" ]]; then
  chkpfx="data."
  chksfx=".dbl"
  goodmsg="> Done"
  tstop_msg="> Last time step - tstop reached"
  maxsteps_msg="> Last time step - maxsteps reached\n\n"
  logfile="pluto.0.log"

elif [[ "$GRID" == "Chombo" ]]; then
  chkpfx="chk."
  chksfx=".3d.hdf5"
# TODO: Check whether this is still correct
  goodmsg="Everything completed"
  tstop_msg="> Last time step - tstop reached"
  maxsteps_msg="> Last time step - maxsteps reached\n\n"
  logfile="pout.0"

else
  echo Unknown grid mode
  exit 1

fi

# Checkpoint file numbering pattern
chknumpat="[0-9]{4}"


#_________________________________________________
# Main bit

# First, change to working directory
cd $WDIR || exit

# Find log and output directories from ini file
OUT_DIR=$(grep -h output_dir $inifile | awk '{print $2}')
OUT_DIR=${OUT_DIR:-"."}
LOG_DIR=$(grep -h log_dir $inifile | awk '{print $2}')
LOG_DIR=${LOG_DIR:-$OUT_DIR}

# Default executable
EXEC=${EXEC:-"./pluto"}

# Default mpirun
MPIRUN=${MPIRUN:-"mpirun"}

# START_NJOB=1 if START_NJOB is undefined because
# it is the first submission
START_NJOB=${START_NJOB:-1}

# NJOB=START_NJOB if NJOB is undefined because
# it is the first submission or a manual restart.
NJOB=${NJOB:-$START_NJOB}

# NCRASH=0 if NCRASH is undefined because
# it is the first submission or a manual restart.
NCRASH=${NCRASH:-0}

# Restart file number used from here on. On initial
# submit or initial manual restart, it is always
# not defined. Upon auto resubmit, it is an empty string
# so the default value here does not take effect.
RESTARTNUM=${RESTARTNUM-$RESTART_FNUM}

# Add inifile to PLUTOARGS
PLUTOARGS+=( -i $inifile )

# If this is a restart
if [[ ($NJOB -gt 1) && ($NJOB -le $NJOBS) ]]; then

  echo "Preparing restart for job $NJOB / $NJOBS"

  # If no initial restartnumber RESTARTNUM is given,
  # grab the number of the last checkpoint file by time
  if [[ -z "$RESTARTNUM" ]]; then

    last_chk_num=$(find $OUT_DIR -name $chkpfx*$chksfx | sort | tail -1 | grep -oE $chknumpat | tail -1)

    # Some checks
    if [[ -z "$last_chk_num" ]]; then
      echo Unable to determine last checkpoint number. Exiting.
      exit 1

    elif [[ "$last_chk_num" == "0000" ]]; then
      echo Cannot use checkpointfile 0000.
      exit 1
    fi

    # This will be the restart number, with leading 0s removed
    RESTARTNUM=$((10#$last_chk_num))

  fi

  # Patch PLUTO arguments to allow restart.
  echo "Check point file number: $RESTARTNUM"
  PLUTOARGS+=( -restart $RESTARTNUM )

fi

# Add maxsteps argument to pluto if MAXSTEPS was given above
if [[ -n "$MAXSTEPS" ]]; then
  PLUTOARGS+=( -maxsteps $MAXSTEPS )
fi

# Add maxtime argument to pluto if MAXTIME was given above
if [[ -n "$MAXTIME" ]]; then
  PLUTOARGS+=( -maxtime $MAXTIME )
fi


# Launch PLUTO
echo "Launching PLUTO for job $NJOB / $NJOBS"
echo "$MPIRUN" "$EXEC" "${PLUTOARGS[@]}"
$MPIRUN $EXEC "${PLUTOARGS[@]}"


# Some housekeeping
# Create a log dir for each submission and copy files into it
logsfx=$(seq -f %02g $NJOB $NJOB)-$PBS_JOBID
cp $LOG_DIR/pluto.0.log $LOG_DIR/pluto.0.$logsfx.log

# Increment crash counter if job didn't finish properly.
if tail $LOG_DIR/$logfile | grep -q "$goodmsg"; then
  echo This run ended successfully.
else
  ((NCRASH++))
  echo This run ended abnormally. Incrementing crash counter.
fi
echo "$NCRASH / $NCRASHMAX crashes have so far occurred."


# If we've done all jobs
if [[ $NJOB -ge $NJOBS ]]; then
  echo "That was the last run. Jobs completed: $NJOB / $NJOBS."
  exit 0

# Exit if more than NCRASHMAX crashes occurred
elif [[ $NCRASH -ge $NCRASHMAX ]]; then
  echo "$NCRASH / $NCRASHMAX abnormal completions occurred. Further restarts cancelled."
  exit 0

# Test whether simulation has finished by either reaching tstop or maxsteps
elif tail $LOG_DIR/$logfile | grep -q "$tstop_msg"; then
  echo Simulation completed: tstop reached.

elif tail $LOG_DIR/$logfile | grep -q "$maxsteps_msg"; then
  echo Simulation completed: maxsteps reached.

# If we still have jobs to do, increment job counter, and resubmit.
else
  ((NJOB++))
  echo "Submitting job $PBS_JOBNAME: number $NJOB / $NJOBS"
  qsub -W depend=afterany:$PBS_JOBID -v NJOB=$NJOB,NCRASH=$NCRASH,RESTARTNUM="" $PBS_JOBNAME

fi


