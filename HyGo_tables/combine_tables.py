import numpy as np

table1 = np.genfromtxt('table_r6.xvg')
table2 = np.genfromtxt('table_r12.xvg')

np.savetxt('tablep.xvg', np.vstack([table1[:,0], table1[:,1], table1[:,2], -table1[:,5], -table1[:,6], table2[:,5], table2[:,6]]).T, fmt='%0.8f')
