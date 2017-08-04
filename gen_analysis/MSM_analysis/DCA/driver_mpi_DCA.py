import os
# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()


script_nm = '/Run_DCA_'+str(rank)+'.sh'
cdir = os.getcwd()
os.system(cdir+script_nm)

