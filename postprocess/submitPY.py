if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 13, 21 ]
	string=open('postprocess.py').readlines() #--- python script
	#---
	nphi = 4
#        MC = [-0.25,-0.1,0.0,0.1,0.2]
#        B = [1.1,1.2,1.3]
#        DF = [1.3,1.4,1.5]
#	PHI = np.logspace(-5.0,-2.0,nphi,endpoint=True)
#	PHI = np.linspace(2.3,2.9,nphi,endpoint=True)
#	PHI = np.linspace(0.05,0.45,nphi,endpoint=True)
	PHI=[0,500000,1000000,1500000,2000000]
	#---
#	nn = 4
#	NTHRESH = np.linspace(0.05,0.11,nn,endpoint=True)
	#---
	jobname = 'd2min'

#	PHI = [[PHI[iphi],NTHRESH[inn]] for iphi in xrange( nphi ) for inn in xrange(nn)]
#	nphi = len(PHI)
	for iphi in xrange( nphi ):
		#---	
		inums = lnums[ 0 ] - 1
		string[ inums ] = "\tjobname  = '%s'\n" % ('%s%s'%(jobname,iphi)) #--- change job name
		#---	densities
		inums = lnums[ 1 ] - 1
#		string[ inums ] = "\targv=\'-p\tmccc\t%e\t-p\tbvall\t%e\t-p\tDfff\t%e\'\n"%(mc, bval, df)
		string[ inums ] = "\targv2nd=\' -p\titime\t%s\'"%(PHI[iphi])

		sfile=open('junk%s.py'%iphi,'w');sfile.writelines(string);sfile.close()
		os.system( 'python junk%s.py'%iphi )
		os.system( 'rm junk%s.py'%iphi )
