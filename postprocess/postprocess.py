def makeOAR( EXEC_DIR, node, core, tpartitionime, PYFIL, argv,argv2nd):
	#--- set environment variables
	sfile = open('.env','w')
	print('%s\n%s'%(argv,argv2nd),file=sfile)
	sfile.close()

	someFile = open( 'oarScript.sh', 'w' )
	print('#!/bin/bash\n',file=someFile)
	print('EXEC_DIR=%s\n'%( EXEC_DIR ),file=someFile)
#	print >> someFile, 'papermill --prepare-only %s/%s ./output.ipynb %s %s'%(EXEC_DIR,PYFIL,argv,argv2nd) #--- write notebook with a list of passed params
	print('jupyter nbconvert --execute $EXEC_DIR/%s --to html --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html'%(PYFIL), file=someFile)
	someFile.close()										  
#
if __name__ == '__main__':
	import os
#
	runs	 = [0,1,2]
	jobname  = {
				'1':'ElasticityT300/CoNiCrFe/itime0', 
				'2':'MlTrain/CoNiCrFeTSNE', 
				}['2']
	DeleteExistingFolder = True
	readPath = os.getcwd() + {
								'1':'/../testRuns/Preparation/ElasticityT300/CoNiCrFe/itime200',
								'2':'/../testRuns/glassCoNiCrFe',
 							}['2'] #--- source
	EXEC_DIR = '.'     #--- path for executable file
	durtn = '23:59:59'
	mem = '128gb'
	partition = ['parallel','cpu2019','bigmem','single'][2] 
	argv = "path=%s"%(readPath) #--- don't change! 
	argv2nd = "itime=2000000\nindx=2" 
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
		os.system( 'cp LammpsPostProcess2nd.py OvitosCna.py %s' % ( writPath ) ) #--- cp python module
		makeOAR( writPath, 1, 1, durtn, PYFIL, argv+"/Run%s"%counter, argv2nd) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh .env %s; cp %s/%s %s' % ( writPath, EXEC_DIR, PYFIL, writPath ) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --mem=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh'\
						   % ( partition, mem, durtn, jobname[:4], counter, jobname[:4], counter, jobname[:4], counter \
						       , writPath, 1, 1, writPath ) ) # --- runs oarScript.sh!
											 

