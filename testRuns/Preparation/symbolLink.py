

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
             'Co5Cr2Fe40Mn27Ni26'
         ] 
   
	for loopi in PHI: #--- composition 
		for loop in np.arange(0,200+1,1): #[50,101,195]: #np.arange(0,200+8,8): #--- itimes
			jobname  = 'ElasticityT300/%s/itime%s'%(loopi,loop)
			job_id = int(open('%s/jobID.txt'%jobname).readlines()[0].split()[3])
			N = 3 #--- nruns
	#---
			job_ids = [ job_id + i for i in xrange( N ) ]
			for id_job, counter in zip( job_ids, xrange( sys.maxint ) ):
				writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
				for file_name in [ 'dump_init.xyz', 'dump_def6.xyz' ]:
					os.system( 'ln -s /scratch/%s/%s %s/' % ( id_job, file_name, writPath ) )
