#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/Preparation/test

MEAM_library_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles

module load mpich/3.2.1-gnu

python pyScript.py

$EXEC_DIR/lmp_serial < in.txt -screen echo -var OUT_PATH . -var MEAM_library_DIR /home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles -var cutoff 2.15443469003
