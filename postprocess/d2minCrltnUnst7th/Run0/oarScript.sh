#!/bin/bash

EXEC_DIR=/Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/d2minCrltnUnst7th/Run0

papermill --prepare-only /Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/d2minCrltnUnst7th/Run0/analyzePlasticity.ipynb ./output.ipynb  -p path '/Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/../BmgData/Run0'  -p itime 2000000
jupyter nbconvert --execute /Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/postprocess/d2minCrltnUnst7th/Run0/output.ipynb --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html
