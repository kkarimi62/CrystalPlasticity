#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/lammps-29Oct20/src
MEAM_library_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles/

module load mpich/3.2.1-gnu
$EXEC_DIR/lmp_serial < in_equilibrate.txt -var pair_coeff_args $MEAM_library_DIR/library_CoNiCrFeMn.meam Co Ni Cr Fe Mn $MEAM_library_DIR/parameters.meam Co Ni Cr Fe Mn
