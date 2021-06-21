#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/Preparation/test7thMelt

MEAM_library_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles

module load mpich/3.2.1-gnu

infile=in.txt
#infile=lattice.in
$EXEC_DIR/lmp_serial < $infile -echo screen -var OUT_PATH . -var MEAM_library_DIR /home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles -var cutoff 1.4
