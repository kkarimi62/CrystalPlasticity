#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/Preparation/test7thMelt

MEAM_library_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles

module load mpich/3.2.1-gnu

$EXEC_DIR/lmp_serial < in.txt -echo screen -var OUT_PATH . -var MEAM_library_DIR /home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles -var cutoff 1.7307086614173226
