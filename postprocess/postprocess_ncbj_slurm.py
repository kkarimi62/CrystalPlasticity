from backports import configparser
def makeOAR( EXEC_DIR, node, core, tpartitionime, PYFIL, argv):
	#--- parse conf. file
	confParser = configparser.ConfigParser()
	confParser.read('config.ini')
	#--- set parameters
	confParser.set('parameters','itime0','0')
	confParser.set('parameters','itime','2000000')
	confParser.set('input files','path',argv)
	confParser.set('input files','fileIndex','5')
	#--- write
	confParser.write(open('config.ini','w'))	
	#--- set environment variables

	someFile = open( 'oarScript.sh', 'w' )
	print('#!/bin/bash\n',file=someFile)
	print('EXEC_DIR=%s\n source /mnt/opt/spack-0.17/share/spack/setup-env.sh\n\nspack load python@3.8.12%%gcc@8.3.0\n\n'%( EXEC_DIR ),file=someFile)
	if convert_to_py:
		print('jupyter nbconvert --to script %s --output py_script\n'%PYFIL,file=someFile)
		print('python3 py_script.py\n',file=someFile)
		 
	else:
		print('jupyter nbconvert --execute $EXEC_DIR/%s --to html --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.allow_errors=True;ls output.html'%(PYFIL), file=someFile)
	someFile.close()										  
#
if __name__ == '__main__':
	import os
#
	runs	 = [0,1,2]
	nNode    = 1
	nThreads = 1
	jobname  = {
				'1':'ElasticityT300/Co5Cr2Fe40Mn27Ni26/eps2/itime0', 
				'2':'MlTrain/Co5Cr2Fe40Mn27Ni26/itime0/Angular', 
				'3':'PairCrltnT300/Co5Cr2Fe40Mn27Ni26', 
				'4':'VorAnlT300/Co5Cr2Fe40Mn27Ni26', 
				'5':'D2minAnalysisT300/Co5Cr2Fe40Mn27Ni26', 
				'6':'MlTrain/granular/itime0', 
				}['1']
	DeleteExistingFolder = True
	readPath = os.getcwd() + {
								'1':'/../testRuns/Preparation/ElasticityT300/Co5Cr2Fe40Mn27Ni26/eps2/itime0',
								'2':'/../testRuns/glassCo5Cr2Fe40Mn27Ni26',
								'3':'/../testRuns/Preparation/CuZrNatom32KT300Tdot1E-1Sheared',
								'4':'/../testRuns/granular/silviaData/DATA_GRAINS/seed1_1001',
 							}['1'] #--- source
	EXEC_DIR = '.'     #--- path for executable file
	durtn = '23:59:59'
	mem = '128gb'
	partition = ['INTEL_PHI'][0] 
	argv = "%s"%(readPath) #--- don't change! 
	PYFILdic = { 
		0:'ElasticConstants.ipynb',
		1:'analyzePlasticity.ipynb',
		3:'junk.ipynb',
		}
	keyno = 1
	convert_to_py = True
#---
#---
	PYFIL = PYFILdic[ keyno ]
	
	#--- update argV
	#---
	if DeleteExistingFolder:
		print('rm %s'%jobname)
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
                                 --chdir %s --ntasks-per-node=%s --nodes=%s %s/oarScript.sh >> jobID.txt'\
                            % ( partition, mem, durtn, jobname.split('/')[0], counter, jobname.split('/')[0], counter, jobname.split('/')[0], counter \
                                , writPath, nThreads, nNode, writPath ) ) # --- runs oarScript.sh! 

