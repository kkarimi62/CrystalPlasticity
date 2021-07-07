def makeOAR( EXEC_DIR, node, core, time, PYFIL ):
	someFile = open( 'oarScript.sh', 'w' )
	print >> someFile, '#!/bin/bash\n'
	print >> someFile, 'EXEC_DIR=%s\n' %( EXEC_DIR )
	print >> someFile, 'MEAM_library_DIR=%s\n' %( MEAM_library_DIR )
	print >> someFile, 'module load mpich/3.2.1-gnu\n'

        '''
	#--- run python script
	pyScript = open( '%s/pyScript.py'%writPath, 'w' )
	print >> pyScript, 'import imp\ngn=imp.load_source(\'generate.name\',\'%s/generate.py\')'%(PYFIL)
	print >> pyScript, 'gn.Generate( %s, %s, %s,title = \'data.txt\',ratio1 = %s, ratio2 = %s, ratio3 = %s, ratio4 = %s, ratio5 = %s )'%(natom, ntypes, rho, 0.05, 0.26, 0.02, 0.4, 0.27)
	pyScript.close()
	print >> someFile, 'python pyScript.py\n'
        '''

	#--- run python script 
	OUT_PATH = '.'
	if SCRATCH:
		OUT_PATH = '/scratch/${SLURM_JOB_ID}'
#	 print >> someFile, "$EXEC_DIR/%s < in.txt -var OUT_PATH %s -var MEAM_library_DIR %s"%( EXEC, OUT_PATH, MEAM_library_DIR )
#	cutoff = 1.0 / rho ** (1.0/3.0)
	if EXEC == 'lmp_serial': 
		print >> someFile, "$EXEC_DIR/%s < in.txt -echo screen -var OUT_PATH %s -var MEAM_library_DIR %s -var cutoff %s -var natom %s"%( EXEC, OUT_PATH, MEAM_library_DIR, cutoff, natom )
	elif EXEC == 'lmp_mpi':
		print >> someFile, "mpirun -np %s $EXEC_DIR/%s < in.txt -echo screen -var OUT_PATH %s -var MEAM_library_DIR %s -var cutoff %s -var natom %s"%( nThreads, EXEC, OUT_PATH, MEAM_library_DIR, cutoff, natom )
	someFile.close()										  


if __name__ == '__main__':
	import os
        import numpy as np

	nruns	 = 1
	nThreads = 32
	jobname  = 'test11thMelt2ndNatom500kQrate0.1-2nd'
#	sourcePath = os.getcwd() + '/dataFiles'
	EXEC_DIR = '/home/kamran.karimi1/Project/git/lammps2nd/lammps/src' #--- path for executable file
	MEAM_library_DIR='/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles' #--- meam potential parameters
	PYFIL = '/home/kamran.karimi1/Project/git/CrystalPlasticity/py'
	EXEC = 'lmp_mpi' #'lmp_serial'
	durtn = '95:59:59'
	SCRATCH = True
	mem = '8gb'
	partition = 'cpu2019' #'cpu2039' #'parallel' #'single' #'parallel'
	#--- sim. parameters
	natom = 50000 * 2 * 2 * 2   
	ntypes = 5
        cutoff = 3.58
    #   cutoffs = np.linspace((1.0-0.5)*cutoff,(1+0.5)*cutoff,nruns)
#	rho = 0.1
	#---
	os.system( 'rm -rf %s' % jobname ) #--- rm existing
	os.system( 'rm jobID.txt' )
	# --- loop for submitting multiple jobs
	counter = 0
	for irun in xrange( nruns ):
#               cutoff = cutoffs[ irun ]
		print ' i = %s' % counter
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		os.system( 'mkdir -p %s' % ( writPath ) ) # --- create folder
		if irun == 0: #--- cp to directory
			path=os.getcwd() + '/%s' % ( jobname)
			os.system( 'cp %s/%s %s' % ( EXEC_DIR, EXEC, path ) ) # --- create folder & mv oar scrip & cp executable
		#---
		os.system( 'cp in_melt.txt %s/in.txt ' % writPath ) #--- lammps script: periodic x, pxx, vy, load
#		os.system( 'cp in_test.txt %s/in.txt ' % writPath ) #--- lammps script: periodic x, pxx, vy, load
		#---
		#---
		makeOAR( path, 1, nThreads, durtn, PYFIL ) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s' % ( writPath) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --mem=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh >> jobID.txt'\
						   % ( partition, mem, durtn, jobname, counter, jobname, counter, jobname, counter \
						       , writPath, nThreads, 1, writPath ) ) # --- runs oarScript.sh! 
		counter += 1
											 
	os.system( 'mv jobID.txt %s' % ( os.getcwd() + '/%s' % ( jobname ) ) )
