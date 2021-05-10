#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/test/Run0

papermill --prepare-only /home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/test/Run0/ElasticConstants.ipynb ./output.ipynb  -p path '/home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/../testRuns/test/Run0' 
jupyter nbconvert --execute /home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/test/Run0/output.ipynb --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html
