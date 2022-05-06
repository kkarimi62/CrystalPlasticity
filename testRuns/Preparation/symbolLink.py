

if __name__ == '__main__':
	import os
	import sys
	import numpy as np
	#--- 
	PHI = [  
#             'FeNi',
#             'CoNiFe',
#             'CoNiCrFe',
#             'CoCrFeMn',
#             'CoNiCrFeMn',
             'Co5Cr2Fe40Mn27Ni26',
#			  'CuZr3'	
         ] 

	EPS = { 
#			'0':1.0,
#			'1':2.0,
#			'2':4.0,
			'3':8.0
		  }
	
   
	runs=[0] #,1,2] #--- nruns
	for loopi in PHI: #--- composition 
		for epsi in EPS: #--- strain 
			for loop in np.arange(0,200+1,8): #[50,101,195]: #np.arange(0,200+8,8): #--- itimes
				jobname  = 'ElasticityT300/%s/eps%s/itime%s'%(loopi,epsi,loop)
#				jobname  = 'ElasticityT300/%s/itime%s'%(loopi,loop)
				job_id = int(open('%s/jobID.txt'%jobname).readlines()[0].split()[3])
	#---
				N = len(runs)
				job_ids = [ job_id + i for i in xrange( N ) ]
				assert len(job_ids) == len(runs)
				for id_job, counter in zip( job_ids, runs): #xrange( sys.maxint ) ):
					writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
					for file_name in ['shearOscillation.xyz']: #[ 'dump_init.xyz', 'dump_def6.xyz' ]:
						os.system( 'ln -s /scratch/%s/%s %s/' % ( id_job, file_name, writPath ) )
