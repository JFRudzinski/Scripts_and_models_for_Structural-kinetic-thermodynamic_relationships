#!/bin/bash

#EXE="/data/isilon/rudzinski/soft_backup/Gromacs_Tables_07.02.13/gmx4_tables-intracut"
EXE="/usr/data/rudzinski/soft_backup2/Gromacs_Tables_07.02.13/gmx4_tables-intracut"

INPUT="force_dh_er-80.0.dat"
OUTPUT="table_NP_NP_er-80.0.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


