#!/bin/bash

#EXE="/data/isilon/rudzinski/soft_backup/Gromacs_Tables_07.02.13/gmx4_tables-intracut"
EXE="/usr/data/rudzinski/soft_backup2/Gromacs_Tables_07.02.13/gmx4_tables-intracut"

INPUT="force_db_rcm-0.5_eps-9_epsssm-0.0495_epsdb-4.95.dat"
OUTPUT="table_CA_CA_eps-9_epsssm-0.0495_epsdb-4.95.xvg"
TYPE="nb"
BASIS="Bspline"
RMAX="20.00"
$EXE $INPUT $OUTPUT $TYPE $BASIS $RMAX


