if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 16, 17 ]
	string=open('testRuns2nd.py').readlines() #--- python script
	#---
	PHI = [	 'FeNi', 
             'CoNiFe', 
             'CoNiCrFe',
             'CoCrFeMn',
             'CoNiCrFeMn',
             'Co5Cr5Fe40Mn27Ni26']

	nphi = len(PHI)
	#---
#	nn = 4
#	NTHRESH = np.linspace(0.05,0.11,nn,endpoint=True)
	#---
	jobname = 'glass'

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
		string[ inums ] = "\targs=\'%s_glass.dump %s_gr.xyz 1\'\n"%(PHI[iphi],PHI[iphi])

		sfile=open('junk%s.py'%iphi,'w');sfile.writelines(string);sfile.close()
		os.system( 'python3 junk%s.py'%iphi )
		os.system( 'rm junk%s.py'%iphi )
