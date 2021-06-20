def makeOAR( EXEC_DIR, node, core, tpartitionime, PYFIL, argv,argv2nd):
	someFile = open( 'oarScript.sh', 'w' )
	print >> someFile, '#!/bin/bash\n'
	print >> someFile, 'EXEC_DIR=%s\n' %( EXEC_DIR )
	print >> someFile, 'papermill --prepare-only %s/%s ./output.ipynb %s %s'%(EXEC_DIR,PYFIL,argv,argv2nd) #--- write notebook with a list of passed params
	print >> someFile, 'jupyter nbconvert --execute %s/output.ipynb --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html'%(EXEC_DIR)
	someFile.close()										  
#
if __name__ == '__main__':
	import os
#
	nruns	 = 1
	jobname  = 'd2min' 
#	readPath = os.getcwd() + '/../testRuns/test8thUnwrapped' # --- source
	readPath = os.getcwd() + '/../BmgData' # --- source
	EXEC_DIR = '.'     #--- path for executable file
	durtn = '23:59:59'
	mem = '8gb' #'128gb'
	partition = 'single' #'cpu2019' #'bigmem' #'single' #'parallel' #'single'
	argv = " -p path \'%s"%(readPath) #--- don't change! 
	argv2nd = " -p itime %s"%(2000000) 
	PYFILdic = { 
		0:'ElasticConstants.ipynb',
		1:'analyzePlasticity.ipynb',
		}
	keyno = 1
#---
#---
	PYFIL = PYFILdic[ keyno ] 
	#--- update argV
	#---
	os.system( 'rm -rf %s' % jobname ) # --- rm existing
	# --- loop for submitting multiple jobs
	counter = 0
	for irun in xrange( nruns ):
		print ' i = %s' % counter
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		os.system( 'mkdir -p %s' % ( writPath ) ) # --- create folder
		os.system( 'cp LammpsPostProcess.py %s' % ( writPath ) ) #--- cp python module
		makeOAR( writPath, 1, 1, durtn, PYFIL, argv+"/Run%s\'"%irun, argv2nd) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s; cp %s/%s %s' % ( writPath, EXEC_DIR, PYFIL, writPath ) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --mem=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh'\
						   % ( partition, mem, durtn, jobname, counter, jobname, counter, jobname, counter \
						       , writPath, 1, 1, writPath ) ) # --- runs oarScript.sh!
		counter += 1
											 

