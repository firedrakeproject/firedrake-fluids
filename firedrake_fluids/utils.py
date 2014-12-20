from mpi4py import MPI

def rank0():
   """ Determine whether or not the calling process has MPI rank 0. """
   if(MPI.COMM_WORLD.Get_rank() == 0):
      return True
   else:
      return False
