#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/Preparation/test

MEAM_library_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles

module load mpich/3.2.1-gnu

python pyScripy.py

$EXEC_DIR/lmp_serial < in.txt -var OUT_PATH . pair_coeff_args $MEAM_library_DIR/library_CoNiCrFeMn.meam Co Ni Cr Fe Mn $MEAM_library_DIR/parameters.meam Co Ni Cr Fe Mn
