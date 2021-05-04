def makeOAR( EXEC_DIR, node, core, time, PYFIL ):
	someFile = open( 'oarScript.sh', 'w' )
	print >> someFile, '#!/bin/bash\n'
	print >> someFile, 'EXEC_DIR=%s\n' %( EXEC_DIR )
	print >> someFile, 'module load mpich/3.2.1-gnu'

	#--- run python script 
#	OUT_PATH = '.'
#	if SCRATCH:
#		OUT_PATH = '/scratch/${SLURM_JOB_ID}'
#	print >> someFile, "mpirun -np %s $EXEC_DIR/%s < in.txt -var boxx %4.3e -var boxz0 %4.3e -var OUT_PATH %s -echo screen -var DATA_FILE \'data.txt\' -var ZHI ${zmax} -var RMAX ${rmax}" %( nThreads, PYFIL, diam, height, OUT_PATH )
	print >> someFile, "$EXEC_DIR/%s < in.txt" %( PYFIL )
	someFile.close()										  


if __name__ == '__main__':
	import os

	nruns	 = 1
	nThreads = 1
	jobname  = 'test'
	sourcePath = os.getcwd() + '/dataFiles'
	EXEC_DIR = '/home/kamran.karimi1/Project/git/CrystalPlasticity/lammps-29Oct20/src' #--- path for executable file
	PYFIL = 'lmp_serial' 
	durtn = '23:59:59'
	SCRATCH = None
	partition = 'single' #'parallel'
	#--- update data.txt and lammps script
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
		os.system( 'cp shearMG300-11.in %s/in.txt ' % writPath ) #--- lammps script: periodic x, pxx, vy, load
		os.system( 'cp %s/FeNi_glass_300.data %s ' % (sourcePath, writPath) ) #--- lammps script: periodic x, pxx, vy, load
		os.system( 'cp %s/parameters.meam %s ' % (sourcePath, writPath) ) 
		os.system( 'cp %s/library_CoNiCrFeMn.meam %s ' % (sourcePath, writPath) ) 
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
