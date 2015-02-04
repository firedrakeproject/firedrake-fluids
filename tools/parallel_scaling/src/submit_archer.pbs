#!/bin/bash --login

#PBS -N __JOBNAME__
#PBS -l walltime=24:00:00
#PBS -l select=__NODES__
#PBS -A n01-IC1
#PBS -M c.jacobs10@imperial.ac.uk

PROJECT=__PROJECT__

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
 
echo Running in $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
 
LOGFILE=${PBS_JOBID}.log
 
module swap PrgEnv-cray PrgEnv-gnu
module unload python
module add anaconda
 
# Add fdrake module path
export FDRAKE_DIR=/work/y07/y07/fdrake
module use $FDRAKE_DIR/modules
 
# Load build environment
module load fdrake-build-env
module load numpy
module load spud
 
export LD_LIBRARY_PATH=${ANACONDA_LIB}:${LD_LIBRARY_PATH}
 
export FIREDRAKE_FFC_KERNEL_CACHE_DIR=$WORK/firedrake-cache
export PYOP2_LAZY=0
export PYOP2_BACKEND_COMPILER=gnu
export PYOP2_CACHE_DIR=$WORK/pyop2-cache
export PYOP2_SIMD_ISA=avx
# Prevent matplotlib from accessing /home
export HOME=$WORK
export XDG_CONFIG_HOME=''

# Paths to local builds
export FIREDRAKE_FLUIDS_PATH=$WORK/firedrake-fluids-bitbucket/
export PYTHONPATH=$WORK/install_new:$WORK/build_new:$WORK/build_new/petsc:$FIREDRAKE_FLUIDS_PATH:$PYTHONPATH
export PYTHONPATH=$WORK/install_new/lib/python2.7/site-packages/:$PYTHONPATH

export PETSC_ARCH=cray-gnu-shared
export PETSC_DIR=$WORK/build_new/petsc

# MPI (man intro_mpi)
export MPICH_NEMESIS_ASYNC_PROGRESS=MC
export MPICH_MAX_THREAD_SAFETY=multiple
export MPICH_CPUMASK_DISPLAY=1
 
# PyOP2 environment variables
export PYOP2_PRINT_SUMMARY=1
export PYOP2_NO_FORK_AVAILABLE=1

echo -n Started at 
date
 
echo
echo Running $PROJECT
echo

#env
#module list

aprun -n __PROCS__ -N __MAXPROCS__ python $FIREDRAKE_FLUIDS_PATH/firedrake_fluids/shallow_water.py $PBS_O_WORKDIR/$PROJECT | tee $LOGFILE

echo -n Finished at 
date
