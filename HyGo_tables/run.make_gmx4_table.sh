#!/bin/bash

EXE="/data/isilon/rudzinski/soft_backup/Gromacs_Tables_07.02.13/gmx4_tables-intracut"

INPUT="force_C_C.dat"
OUTPUT="table_C_C.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="force_N_N.dat"
OUTPUT="table_N_N.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="force_N_C.dat"
OUTPUT="table_N_C.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="force_O_O.dat"
OUTPUT="table_O_O.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="force_N_O.dat"
OUTPUT="table_N_O.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX

INPUT="force_C_O.dat"
OUTPUT="table_C_O.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


