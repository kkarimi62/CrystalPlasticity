

if __name__ == '__main__':
	import os
	import sys
	#--- 
        for loop in range(0,200,1):
            jobname  = 'Co5Cr2Fe40Mn27Ni26T300Elasticity%s'%loop
            job_id = int(open('%s/jobID.txt'%jobname).readlines()[0].split()[3])
#            job_id = 12158618 
            N = 3
	#---
	    job_ids = [ job_id + i for i in xrange( N ) ]
	    for id_job, counter in zip( job_ids, xrange( sys.maxint ) ):
		writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
		for file_name in [ 'dump_init.xyz', 'dump_def6.xyz' ]:
			os.system( 'ln -s /scratch/%s/%s %s/' % ( id_job, file_name, writPath ) )
