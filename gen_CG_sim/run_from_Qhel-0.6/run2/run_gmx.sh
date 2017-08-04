#!/bin/bash

#####################################################################

#System Specific
sys="ala7"
######################################################################

#Directories
bin="/sw/linux/gromacs/4.5/full/4.5.3/bin"

######################################################################

#Executables
grompp="$bin/grompp"
mdrun="$bin/mdrun"
tpbconv="$bin/tpbconv"

######################################################################

#File Names
mdp="${sysT}.mdp"
mdout="mdout-${sys}.mdp"
gro="../${sys}_Qhel-0.6.gro"
top="${sys}.top"
edr="${sysT}.edr"
tpr="${sysT}.tpr"
oldtrr="${sysT}.PAR.trr"
trr="${sysT}.trr"
xtc="${sysT}.xtc"
ndx="index_types.ndx"
gro_out="${sysT}.confout.gro"
md_log="mdlog-${sysT}.log"
mpirun="mpirun"
gr_log="grompp.log"
md_log="mdrun.log"
cpi="state_old.cpt"
time="1000"
frame="4641"
table="table.xvg"

#####################################################################

# Determine processors for mpirun
NODELIST=$( cat $PBS_NODEFILE )
echo Processors used
count=0
for i in $NODELIST
do
    let count=$count+1
    echo Processor: $count $i
done
nnodes=$count

nnodes=8
#####################################################################

#grompp
#Normal
$grompp -f $sys -c $gro -p $sys -n $ndx -o $sys -po $mdout\
        >& $gr_log

#Extend
#$grompp -f $mdp -c $oldtpr -p $top -o $newtpr -po $mdout\
#        >& $gr_log

#$tpbconv -s $newtpr -extend $time -o $new2tpr 

#use trr file for starting conf and velocities
#$grompp -f $mdp -c $gro -t $oldtrr -n index.ndx -time $frame -p $top -o $tpr -po $mdout\
#        >& $gr_log

#Tables
#$grompp -f $mdp -c $gro -p $top -o $tpr -po $mdout\
#        -n $ndx >& $gr_log

#Suffle
#$grompp -f $mdp -c $gro -p $top -o $tpr -po $mdout\
#        -shuffle  -deshuf $gr_des >& $gr_log

#####################################################################

#mdrun
#normal
#$MPIEXEC 
#$mdrun -s $tpr -o $trr -x $xtc -c $gro_out -e $edr -g $md_log >& $md_log

$mdrun -s $sys -o $sys -x $sys -c $gro_out -e $sys -g $md_log >& $md_log

#normal-old-CGmpirun
#$mpirun $mdrun -pd -s $tpr -o $trr -c $gro_out -table $table -tableb $table -n $nnodes -e $edr -g $md_log >& $md_log --mca pml ob1

#extend
#$mpirun $mdrun -s $new2tpr -cpi $cpi -o $trr -c $gro_out -n $nnodes -e $edr -g $md_log >& $md_log





