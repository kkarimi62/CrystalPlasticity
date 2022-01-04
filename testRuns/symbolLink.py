

if __name__ == '__main__':
	import os
	import sys
	#--- 
	jobname  = 'test3rd'
	job_id = 9489151 
	N = 1
	#---
	job_ids = [ job_id + i for i in xrange( N ) ]
	for id_job, counter in zip( job_ids, xrange( sys.maxint ) ):
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		for file_name in [ 'virialStress.xyz' ]:
			os.system( 'ln -s /scratch/%s/%s %s/' % ( id_job, file_name, writPath ) )
