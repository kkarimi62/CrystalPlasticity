#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/lammps-29Oct20/src

module load mpich/3.2.1-gnu
$EXEC_DIR/lmp_serial < in_equilibrate.txt -var pair_coeff_args 'library_CoNiCrFeMn.meam Co Ni Cr Fe Mn parameters.meam Co Ni Cr Fe Mn'
