#!/bin/bash

EXE="/data/isilon/rudzinski/soft_backup/Gromacs_Tables_07.02.13/gmx4_tables-intracut"

INPUT="force_r6.dat"
OUTPUT="table_r6.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="force_r12.dat"
OUTPUT="table_r12.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX
