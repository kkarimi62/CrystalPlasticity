if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 18, 20 ]
	string=open('postprocess.py').readlines() #--- python script
	#---
     PHI = [  'FeNi',
              'CoNiFe',
              'CoNiCrFe',
              'CoCrFeMn',
              'CoNiCrFeMn',
              'Co5Cr5Fe40Mn27Ni26']

#	PHI=range(0,2000000+10000,4*10000) 
	nphi = len(PHI)
	#---
#	nn = 4
#	NTHRESH = np.linspace(0.05,0.11,nn,endpoint=True)
	#---
	jobname = 'Cna'

#	PHI = [[PHI[iphi],NTHRESH[inn]] for iphi in xrange( nphi ) for inn in xrange(nn)]
#	nphi = len(PHI)
	for iphi in range( nphi ):
		#---	
		inums = lnums[ 0 ] - 1
		string[ inums ] = "\tjobname  = '%s'\n" % ('%s%s'%(jobname,PHI[iphi])) #--- change job name
		#---	densities
		inums = lnums[ 1 ] - 1
#		string[ inums ] = "\targv=\'-p\tmccc\t%e\t-p\tbvall\t%e\t-p\tDfff\t%e\'\n"%(mc, bval, df)
#		string[ inums ] = "\targv2nd=\' -p\titime\t%s\'\n"%(PHI[iphi])
		string[ inums ] = "\treadPath = os.getcwd() + \'/../testRuns/glass%s\'\n"%(PHI[iphi])

		sfile=open('junk%s.py'%iphi,'w');sfile.writelines(string);sfile.close()
		os.system( 'python3 junk%s.py'%iphi )
		os.system( 'rm junk%s.py'%iphi )
