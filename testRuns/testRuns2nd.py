def makeOAR( EXEC_DIR, node, core, time):
	someFile = open( 'oarScript.sh', 'w' )
	someFile.write('#!/bin/bash\n')
	someFile.write('EXEC_DIR=%s\n' %( EXEC_DIR ))
#	someFile.write('module unload gcc/4.6.3;module load gcc/9.3.0\n')

	someFile.write("ovitos $EXEC_DIR/OvitosCna.py %s" %( args )) #--- cna analysis in ovito!
	someFile.close()										  


if __name__ == '__main__':
	import os

	nruns	 = 3
	nThreads = 1
	jobname  = 'glassFeNi' #--- existing directory
	args = 'FeNi_glass.dump FeNi_gr.xyz 1' #--- input, output
	sourcePath = os.getcwd() + '/dataFiles'
	EXEC_DIR = '/mnt/home/kkarimi/Project/git/CrystalPlasticity/py' #--- path for executable file
	durtn = '00:59:59'
	SCRATCH = None
	partition = 'single' #'parallel'
	#--- update data.txt and lammps script
	#---
	os.system( 'rm jobID.txt' )
	# --- loop for submitting multiple jobs
	counter = 0
	for irun in range( nruns ):
		print(' i = %s' % counter)
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		if irun == 0: #--- cp to directory
			path=os.getcwd() + '/%s' % ( jobname)
		#---
		#---
		makeOAR( EXEC_DIR, 1, nThreads, durtn) # --- make oar script
		os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s' % ( writPath) ) # --- create folder & mv oar scrip & cp executable
		os.system( 'sbatch --partition=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
						    --chdir %s -c %s -n %s %s/oarScript.sh >> jobID.txt'\
						   % ( partition, durtn, jobname, counter, jobname, counter, jobname, counter \
						       , writPath, nThreads, 1, writPath ) ) # --- runs oarScript.sh! 
		counter += 1
											 
	os.system( 'mv jobID.txt %s' % ( os.getcwd() + '/%s' % ( jobname ) ) )
