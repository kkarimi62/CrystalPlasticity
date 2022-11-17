from backports import configparser
def makeOAR( EXEC_DIR, node, core, tpartitionime, PYFIL, argv):
	#--- parse conf. file
	confParser = configparser.ConfigParser()
	confParser.read('config.ini')
	#--- set parameters
	confParser.set('parameters','itime0','5')
	confParser.set('parameters','itime','15')
	confParser.set('input files','path',argv)
	#--- write
	confParser.write(open('configuration.ini','w'))	
	#--- set environment variables

	someFile = open( 'oarScript.sh', 'w' )
	print('#!/bin/bash\n',file=someFile)
	print('EXEC_DIR=%s\n module load python/anaconda3-2019.10-tensorflowgpu'%( EXEC_DIR ),file=someFile)
#	print >> someFile, 'papermill --prepare-only %s/%s ./output.ipynb %s %s'%(EXEC_DIR,PYFIL,argv,argv2nd) #--- write notebook with a list of passed params
	print('jupyter nbconvert --execute $EXEC_DIR/%s --to html --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html'%(PYFIL), file=someFile)
	someFile.close()										  
#
if __name__ == '__main__':
	import os
#
	runs	 = [0]  #,1,2]
	jobname  = {
				'1':'ElasticityT300/Co5Cr2Fe40Mn27Ni26/itime0', 
				'2':'MlTrain/Co5Cr2Fe40Mn27Ni26/itime0/Angular', 
				'3':'PairCrltnT300/Co5Cr2Fe40Mn27Ni26', 
				'4':'VorAnlT300/Co5Cr2Fe40Mn27Ni26', 
				'5':'D2minAnalysisT300/Co5Cr2Fe40Mn27Ni26', 
				'6':'MlTrain/granular/itime0', 
				}['6']
	DeleteExistingFolder = False
	readPath = os.getcwd() + {
								'1':'/../testRuns/Preparation/ElasticityT300/Co5Cr2Fe40Mn27Ni26/itime0',
								'2':'/../testRuns/glassCo5Cr2Fe40Mn27Ni26',
								'3':'/../testRuns/Preparation/CuZrNatom32KT300Tdot1E-1Sheared',
								'4':'/../testRuns/granular/silviaData/DATA_GRAINS/seed1_1001',
 							}['4'] #--- source
	EXEC_DIR = '.'     #--- path for executable file
	durtn = '23:59:59'
	mem = '128gb'
	partition = ['parallel','cpu2019','bigmem','single'][3] 
	argv = "path=%s"%(readPath) #--- don't change! 
	PYFILdic = { 
		0:'ElasticConstants.ipynb',
		1:'analyzePlasticity.ipynb',
		2:'test2nd.ipynb',
		3:'junk.ipynb',
		}
	keyno = 1
#---
#---
	PYFIL = PYFILdic[ keyno ] 
	#--- update argV
	#---
	if DeleteExistingFolder:
		os.system( 'rm -rf %s' % jobname ) # --- rm existing
	# --- loop for submitting multiple jobs
	for counter in runs:
		print(' i = %s' % counter)
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		os.system( 'mkdir -p %s' % ( writPath ) ) # --- create folder
		os.system( 'cp utility.py LammpsPostProcess2nd.py OvitosCna.py %s' % ( writPath ) ) #--- cp python module
		makeOAR( writPath, 1, 1, durtn, PYFIL, argv+"/Run%s"%counter) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s; cp config.ini %s;cp %s/%s %s' % ( writPath, writPath, EXEC_DIR, PYFIL, writPath ) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --mem=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh'\
						   % ( partition, mem, durtn, jobname[:4], counter, jobname[:4], counter, jobname[:4], counter \
						       , writPath, 1, 1, writPath ) ) # --- runs oarScript.sh!
											 

