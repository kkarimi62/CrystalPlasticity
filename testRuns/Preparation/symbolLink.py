

if __name__ == '__main__':
	import os
	import sys
	#--- 
	jobname  = 'CoNiCrFeT300Elasticity0'
	job_id = 12158618 
	N = 3
	#---
	job_ids = [ job_id + i for i in xrange( N ) ]
	for id_job, counter in zip( job_ids, xrange( sys.maxint ) ):
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		for file_name in [ 'dump_init.xyz', 'dump_def6.xyz' ]:
			os.system( 'ln -s /scratch/%s/%s %s/' % ( id_job, file_name, writPath ) )
