#!/bin/bash

EXEC_DIR=/home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/test3rd1/Run0

papermill --prepare-only /home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/test3rd1/Run0/ElasticConstants.ipynb ./output.ipynb  -p path '/home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/../testRuns/test3rd/Run0' 
jupyter nbconvert --execute /home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess/test3rd1/Run0/output.ipynb --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html
