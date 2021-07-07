if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 35, 36 ]
	string=open('Preparation.py').readlines() #--- python script
	#---
#        MC = [-0.25,-0.1,0.0,0.1,0.2]
#        B = [1.1,1.2,1.3]
#        DF = [1.3,1.4,1.5]
#	PHI = np.logspace(-5.0,-2.0,nphi,endpoint=True)
#	PHI = np.linspace(2.3,2.9,nphi,endpoint=True)
#	PHI = np.linspace(0.05,0.45,nphi,endpoint=True)
	PHI=range(2,42,2) #[1,2,4,6,8,10,12,14,16,18,20,22,24,26,28] #range(0,2000000+10000,4*10000) 
#	PHI = [0, 500000,1000000,1500000,2000000]
	nphi = len(PHI)
	#---
#	nn = 4
#	NTHRESH = np.linspace(0.05,0.11,nn,endpoint=True)
	#---
	jobname = 'testCpuRunTime'

#	PHI = [[PHI[iphi],NTHRESH[inn]] for iphi in xrange( nphi ) for inn in xrange(nn)]
#	nphi = len(PHI)
	for iphi in range( nphi ):
		#---	
		inums = lnums[ 0 ] - 1
		string[ inums ] = "\tnThreads=%s\n"%(PHI[iphi])
		#---	densities
		inums = lnums[ 1 ] - 1
		string[ inums ] = "\tjobname  = '%s'\n" % ('%s%s'%(jobname,iphi)) #--- change job name
#		string[ inums ] = "\targv=\'-p\tmccc\t%e\t-p\tbvall\t%e\t-p\tDfff\t%e\'\n"%(mc, bval, df)
#		string[ inums ] = "\targv2nd=\' -p\titime\t%s\'\n"%(PHI[iphi])
#		string[ inums ] = "\targv2nd=\'itime=%s\'\n"%(PHI[iphi])

		sfile=open('junk%s.py'%iphi,'w');sfile.writelines(string);sfile.close()
		os.system( 'python junk%s.py'%iphi )
		os.system( 'rm junk%s.py'%iphi )
