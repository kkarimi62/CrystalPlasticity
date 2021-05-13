def makeOAR( EXEC_DIR, node, core, time, PYFIL ):
	someFile = open( 'oarScript.sh', 'w' )
	print >> someFile, '#!/bin/bash\n'
	print >> someFile, 'EXEC_DIR=%s\n' %( EXEC_DIR )
	print >> someFile, 'module load mpich/3.2.1-gnu'

	#--- run python scri
	print >> open( 'pyScripy.py', 'w' ), 'import generate as gn; \ 
				   						Generate( %s, %s,(%s,%s), (%s, %s), (%s, %s), \
			 							title = \'data.txt\', \
			 							ratio1 = %s, ratio2 = %s, ratio3 = %s, ratio4 = %s, ratio5 = %s )' \
										%(natom, ntypes, xlo, xhi, ylo, yhi, zlo, zhi, 0.05, 0.26, 0.02, 0.4, 0.27)
	print >> someFile, 'python pyScripy.py'

	#--- run python script 
	OUT_PATH = '.'
	if SCRATCH:
		OUT_PATH = '/scratch/${SLURM_JOB_ID}'
	print >> someFile, "$EXEC_DIR/%s < in.txt -var OUT_PATH %s pair_coeff_args %s/library_CoNiCrFeMn.meam Co Ni Cr Fe Mn %s/parameters.meam Co Ni Cr Fe Mn" %( PYFIL, OUT_PATH, MEAM_library_DIR, MEAM_library_DIR )
	someFile.close()										  


if __name__ == '__main__':
	import os

	nruns	 = 1
	nThreads = 1
	jobname  = 'test'
#	sourcePath = os.getcwd() + '/dataFiles'
	EXEC_DIR = '/home/kamran.karimi1/Project/git/CrystalPlasticity/lammps-29Oct20/src' #--- path for executable file
	MEAM_library_DIR='/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles' #--- meam potential parameters
	PYFIL = 'lmp_serial' 
	durtn = '47:59:59'
	SCRATCH = None
	partition = 'single' #'parallel'
	#--- sim. parameters
	(xlo, xhi ) = (0.0, 1.0)
	(ylo, yhi ) = (0.0, 1.0)
	(zlo, zhi ) = (0.0, 1.0)
	natom = 1000
	ntypes = 5
	#---
	os.system( 'rm -rf %s' % jobname ) #--- rm existing
	os.system( 'rm jobID.txt' )
	# --- loop for submitting multiple jobs
	counter = 0
	for irun in xrange( nruns ):
		print ' i = %s' % counter
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		os.system( 'mkdir -p %s' % ( writPath ) ) # --- create folder
		if irun == 0: #--- cp to directory
			path=os.getcwd() + '/%s' % ( jobname)
			os.system( 'cp %s/%s %s' % ( EXEC_DIR, PYFIL, path ) ) # --- create folder & mv oar scrip & cp executable
		#---
		os.system( 'cp in_equilibrate.txt %s/in.txt ' % writPath ) #--- lammps script: periodic x, pxx, vy, load
		#---
		#---
		makeOAR( path, 1, nThreads, durtn, PYFIL ) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s' % ( writPath) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh >> jobID.txt'\
						   % ( partition, durtn, jobname, counter, jobname, counter, jobname, counter \
						       , writPath, nThreads, 1, writPath ) ) # --- runs oarScript.sh! 
		counter += 1
											 
	os.system( 'mv jobID.txt %s' % ( os.getcwd() + '/%s' % ( jobname ) ) )
