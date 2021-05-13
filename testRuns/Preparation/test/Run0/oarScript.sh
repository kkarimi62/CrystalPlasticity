#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/Preparation/test

MEAM_library_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles

module load mpich/3.2.1-gnu

python pyScript.py

$EXEC_DIR/lmp_serial < in.txt -var OUT_PATH . -var pairCoeffArgs /home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles/library_CoNiCrFeMn.meam Co Ni Cr Fe Mn /home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles/parameters.meam Co Ni Cr Fe Mn
