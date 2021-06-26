#!/bin/bash

EXEC_DIR=/Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/d2minFeNi/Run0

conda activate test-env;papermill --prepare-only /Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/d2minFeNi/Run0/analyzePlasticity.ipynb ./output.ipynb  -p path '/Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/../testRuns/glassFeNi/Run0' 
jupyter nbconvert --execute /Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/d2minFeNi/Run0/output.ipynb --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html
