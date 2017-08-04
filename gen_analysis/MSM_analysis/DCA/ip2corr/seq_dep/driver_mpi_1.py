import os
# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()


script_nm = ['MSM_from_T3res_iter-solve_recursloops_memsave.py', 'Get_HMM_AB-Ass0.py', 'Get_HMM_AB-Bss0.py', 'Get_HMM_BBss0.py' ]
for i in range(size):
    if (i == rank and i == 0 ):
        os.system('/home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/python '+script_nm[i])

