

if __name__ == '__main__':
	import os
	import sys
	#--- 
	jobname  = 'uniax17'
	job_id = 9358490
	N = 3
	#---
	job_ids = [ job_id + i for i in xrange( N ) ]
	for id_job, counter in zip( job_ids, xrange( sys.maxint ) ):
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		for file_name in [ 'cap_stress.dump', 'pair_stress.dump' ]:
			os.system( 'ln -s /scratch/%s/%s %s/' % ( id_job, file_name, writPath ) )
