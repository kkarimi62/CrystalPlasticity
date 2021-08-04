if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 18, 20,26 ]
	string=open('postprocess.py').readlines() #--- python script
	#---
	PHI = [	 'FeNi', 
             'CoNiFe', 
             'CoNiCrFe',
             'CoCrFeMn',
             'CoNiCrFeMn',
             'Co5Cr5Fe40Mn27Ni26']
#	PHI = [1800,1933,2066,2200,2333]
#	PHI=range(0,2000000+10000,4*100000) 
	nphi = len(PHI)
	#---
#	nn = 4
#	NTHRESH = np.linspace(0.05,0.11,nn,endpoint=True)
	#---
	jobname = ''

#	PHI = [[PHI[iphi],NTHRESH[inn]] for iphi in xrange( nphi ) for inn in xrange(nn)]
#	nphi = len(PHI)
	for iphi in range( nphi ):
		#---	
		inums = lnums[ 0 ] - 1
#		string[ inums ] = "\tjobname=\'Natom250kQrate0.5Tstop%sSrate-5Strs'\n" % (PHI[iphi]) #--- change job name
		string[ inums ] = "\tjobname=\'Cna%s'\n" % (PHI[iphi]) #--- change job name
		#---	densities
		inums = lnums[ 1 ] - 1
#		string[ inums ] = "\treadPath = os.getcwd() + \'/../testRuns/Preparation/Natom250kQrate0.5Tstop%sSrate-5\'\n"%(PHI[iphi])
		string[ inums ] = "\treadPath = os.getcwd() + \'/../testRuns/glass%s\'\n"%(PHI[iphi])
#
		inums = lnums[ 2 ] - 1
		string[ inums ] = "\targv2nd = \'indx=%s\'\n"%(iphi)

		sfile=open('junk%s.py'%iphi,'w');sfile.writelines(string);sfile.close()
		os.system( 'python3 junk%s.py'%iphi )
		os.system( 'rm junk%s.py'%iphi )
