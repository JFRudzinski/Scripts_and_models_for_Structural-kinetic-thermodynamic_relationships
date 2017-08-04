python='/home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/python'

cd ip2corr/seq_dep/
${python} get_T3res_CG_seqdep.py
${python} MSM_from_T3res_iter-solve_recursloops_memsave.py
cd ../../

