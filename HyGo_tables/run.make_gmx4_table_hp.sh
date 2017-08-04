#!/bin/bash

#EXE="/data/isilon/rudzinski/soft_backup/Gromacs_Tables_07.02.13/gmx4_tables-intracut"
EXE="/usr/data/rudzinski/soft_backup2/Gromacs_Tables_07.02.13/gmx4_tables-intracut"

INPUT="force_hp_eps-6.03.dat"
OUTPUT="table_CB_CB_eps-6.03_eq_0.67times9epsnc.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


