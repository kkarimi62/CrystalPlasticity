def makeOAR( EXEC_DIR, node, core, time, PYFIL ):
	someFile = open( 'oarScript.sh', 'w' )
	print >> someFile, '#!/bin/bash\n'
	print >> someFile, 'EXEC_DIR=%s\n' %( EXEC_DIR )
	print >> someFile, 'MEAM_library_DIR=%s\n' %( MEAM_library_DIR )
	print >> someFile, 'module load mpich/3.2.1-gnu'

	#--- run python script 
	pyScript = open( '%s/pyScript.py'%writPath, 'w' )
	print >> pyScript, 'import imp\ngn=imp.load_source(\'generate.name\',\'%s/generate.py\')'%(PYFIL)
	print >> pyScript, 'gn.GenerateDataFromDump( \'%s/Run%s/%s\',ntype=5,title = \'data.txt\')'%(sourcePath, irun, 'AsQuenched.dump') #--- generate data file from a dump file
	pyScript.close()
	print >> someFile, 'python pyScript.py\n'

	#--- parse sample height
#	print >> someFile, "IFS=\',\' read -ra a <<< \"$string\"" #--- split into array
#	print >> someFile, "zmax=${a[0]} #--- height"
#	print >> someFile, "rmax=${a[1]} #--- max bond length"
#	print >> someFile, "echo ${zmax}"
#	print >> someFile, "echo${rmax}"

	#--- run lammps: compression test
	OUT_PATH = '.'
	if SCRATCH:
		OUT_PATH = '/scratch/${SLURM_JOB_ID}'
#	print >> someFile, "mpirun -np %s $EXEC_DIR/%s < in.txt -var OUT_PATH %s -var MEAM_library_DIR %s -echo screen" %( nThreads, EXEC, OUT_PATH, MEAM_library_DIR )
	if EXEC == 'lmp_serial': 
		print >> someFile, "$EXEC_DIR/%s < in.txt -var OUT_PATH %s -var MEAM_library_DIR %s -echo screen" %( EXEC, OUT_PATH, MEAM_library_DIR )
	elif EXEC == 'lmp_mpi':
		print >> someFile, "mpirun -np %s $EXEC_DIR/%s < in.txt -var OUT_PATH %s -var MEAM_library_DIR %s"%( nThreads, EXEC, OUT_PATH, MEAM_library_DIR)
	someFile.close()										  


if __name__ == '__main__':
	import os

	nruns	 = 1
	nThreads = 32
	jobname  = 'test11thMelt2ndNatom50kQrate1Srate6'
	sourcePath = os.getcwd() + '/test11thMelt2ndNatom50kQrate1' #--- must be different than sourcePath
	EXEC_DIR = '/home/kamran.karimi1/Project/git/lammps2nd/lammps/src' #--- path for executable file
	MEAM_library_DIR='/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/dataFiles' #--- meam potential parameters
	PYFIL = '/home/kamran.karimi1/Project/git/CrystalPlasticity/py'
	EXEC = 'lmp_mpi' #'lmp_serial' 
	durtn = '167:59:59'
	SCRATCH = None
	partition = 'cpu2019' #'parallel' #'single' #'bigmem' #'parallel'
	mem = '8gb'
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
			os.system( 'cp %s/%s %s' % ( EXEC_DIR, EXEC, path ) ) # --- create folder & mv oar scrip & cp executable
		#---
		os.system( 'cp in_shear.txt %s/in.txt ' % writPath ) #--- lammps script: periodic x, pxx, vy, load
#		os.system( 'cp data.txt %s' % writPath ) #--- py script add vels 
		#---
		makeOAR( path, 1, nThreads, durtn, PYFIL ) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s' % ( writPath) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --mem=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh >> jobID.txt'\
						   % ( partition, mem, durtn, jobname, counter, jobname, counter, jobname, counter \
						       , writPath, nThreads, 1, writPath ) ) # --- runs oarScript.sh! 
		counter += 1
											 
	os.system( 'mv jobID.txt %s' % ( os.getcwd() + '/%s' % ( jobname ) ) )
