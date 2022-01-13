if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 23, 28,35 ]
	string=open('postprocess.py').readlines() #--- python script
	#---
#	PHI=range(0,200,1)
	PHI = [	 
            'FeNi',
            'CoNiFe',
            'CoNiCrFe',
            'CoCrFeMn',
            'CoNiCrFeMn',
            'Co5Cr2Fe40Mn27Ni26'
		]
	nphi = len(PHI)
	#---
	itime = 0	

#	PHI = [[PHI[iphi],NTHRESH[inn]] for iphi in xrange( nphi ) for inn in xrange(nn)]
#	nphi = len(PHI)
	for iphi in range( nphi ):
		#---	
		inums = lnums[ 0 ] - 1
#		string[ inums ] = "\t\'1\':\'ElasticityT300/%s/itime%s',\n" % (PHI[iphi],itime) #--- change job name
#		string[ inums ] = "\t\'3\':\'PairCrltnT300/%s/itime%s',\n" % (PHI[iphi],itime) #--- change job name
		string[ inums ] = "\t\'4\':\'VorAnlT300/%s',\n" % (PHI[iphi]) #--- change job name
		#---	densities
		inums = lnums[ 1 ] - 1
#		string[ inums ] = "\t\'1\':\'/../testRuns/Preparation/ElasticityT300/%s/itime%s\',\n"%(PHI[iphi],itime)
		string[ inums ] = "\t\'2\':\'/../testRuns/glass%s\',\n"%(PHI[iphi])
#
		inums = lnums[ 2 ] - 1
		string[ inums ] = "\targv2nd = \'itime=%s\\nindx=%s\'\n"%(itime*10000,iphi) #(PHI[iphi]*10000)

		sfile=open('junk%s.py'%iphi,'w');sfile.writelines(string);sfile.close()
		os.system( 'python3 junk%s.py'%iphi )
		os.system( 'rm junk%s.py'%iphi )
