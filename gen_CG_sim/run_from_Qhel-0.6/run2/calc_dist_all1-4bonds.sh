#!/bin/bash

#####################################################################

#System Specific
sys="ala7"

blen=0.6
tol=5.0

######################################################################

#Directories
bin="/gpfs/work/jfr148/pkg/gromacs-4.5.3-dp/bin"
suff="-xj_sn-4"

######################################################################

#Executables
g_rdf="${bin}/g_rdf${suff}"
g_bond="${bin}/g_bond${suff}"
g_angle="${bin}/g_angle${suff}"
g_rmsd="${bin}/g_rmsdist${suff}"
g_rmsf="${bin}/g_rmsf${suff}"

#g_rdf="g_rdf4.5.3"
#g_bond="g_bond4.5.3"
#g_angle="g_angle4.5.3"
#g_rmsd="g_rmsdist4.5.3"
#g_rmsf="g_rmsf4.5.3"

######################################################################

#File Names
tpr="${sys}.tpr"
trr="${sys}.xtc"
ndx="index_dist_all1-4bonds.ndx"
native="${sys}.gro"

######################################################################

#Calculate Distributions

#bonds
frames=50001

g_bond -f ${trr} -s ${tpr} -blen ${blen} -tol ${tol} -n ${ndx} -d distance_bonds_1-4.xvg -o bonds_1-4.xvg <<-EOF
       0
EOF

tail -n ${frames} distance_bonds_1-4.xvg &> tmp.xvg
mv tmp.xvg distance_bonds_1-4.xvg

g_bond -f ${trr} -s ${tpr} -blen ${blen} -tol ${tol} -n ${ndx} -d distance_bonds_2-5.xvg -o bonds_2-5.xvg <<-EOF
       1
EOF

tail -n ${frames} distance_bonds_2-5.xvg &> tmp.xvg
mv tmp.xvg distance_bonds_2-5.xvg

g_bond -f ${trr} -s ${tpr} -blen ${blen} -tol ${tol} -n ${ndx} -d distance_bonds_3-6.xvg -o bonds_3-6.xvg <<-EOF
       2
EOF

tail -n ${frames} distance_bonds_3-6.xvg &> tmp.xvg
mv tmp.xvg distance_bonds_3-6.xvg

g_bond -f ${trr} -s ${tpr} -blen ${blen} -tol ${tol} -n ${ndx} -d distance_bonds_4-7.xvg -o bonds_4-7.xvg <<-EOF
       3
EOF

tail -n ${frames} distance_bonds_4-7.xvg &> tmp.xvg
mv tmp.xvg distance_bonds_4-7.xvg

g_gyrate -f ${trr} -s ${tpr} -n index_types.ndx -o Rg.xvg <<-EOF
       0
EOF

tail -n ${frames} Rg.xvg &> tmp.xvg
mv tmp.xvg Rg.xvg

g_rama -f ${trr} -s ../../ala7.for-rama.tpr -o rama.xvg

grep 'ALA-2' rama.xvg &> rama_ALA2.xvg
grep 'ALA-3' rama.xvg &> rama_ALA3.xvg
grep 'ALA-4' rama.xvg &> rama_ALA4.xvg
grep 'ALA-5' rama.xvg &> rama_ALA5.xvg
grep 'ALA-6' rama.xvg &> rama_ALA6.xvg
