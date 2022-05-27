if __name__ == '__main__':
	import sys
	import os
	import numpy as np
	#---
	lnums = [ 24, 37 ]
	string=open('postprocess.py').readlines() #--- python script
	#---
	PHI ={ 
#            '0':'FeNi',
#            '1':'CoNiFe',
#           '2':'CoNiCrFe',
#           '3' :'CoCrFeMn',
#            '4':'CoNiCrFeMn',
#            '5':'Co5Cr2Fe40Mn27Ni26'
			'6':'cuzr',
		}


	EPS = {
#             '0':1.0e-03,
#             '1':0.5e-03,
             '2':0.25e-03,
#			'3':2.0e-02,
#              '4':16.0,
#              '5':32.0,
#              '6':64.0
          }


	nphi = len(PHI)
	#---
	times=np.arange(0,100+1)  #--- run GetFrames(lmpData,times=lmpData.coord_atoms_broken.keys())
	nevery = 10000 #--- run list(lmpData.coord_atoms_broken.keys())[1]

#	PHI = [[PHI[iphi],NTHRESH[inn]] for iphi in xrange( nphi ) for inn in xrange(nn)]
#	nphi = len(PHI)
	count = 0
	for key in PHI:
		for epsi in EPS:
			for itime in times:
			#---	
				inums = lnums[ 0 ] - 1
				string[ inums ] = "\t\'6\':\'MlTrain/%s/itime%s',\n" % (PHI[key],itime) #--- change job name
#				string[ inums ] = "\t\'1\':\'ElasticityT300/%s/eps%s/itime%s',\n" % (PHI[key],epsi,itime) #--- change job name
#			string[ inums ] = "\t\'3\':\'PairCrltnT300/%s',\n" % (PHI[iphi]) #--- change job name
	#		string[ inums ] = "\t\'4\':\'VorAnlT300/%s',\n" % (PHI[iphi]) #--- change job name
	#		string[ inums ] = "\t\'5\':\'D2minAnalysisT300/%s',\n" % (PHI[iphi]) #--- change job name
		#---	densities
				inums = lnums[ 1 ] - 1
#				string[ inums ] = "\t\'1\':\'/../testRuns/Preparation/ElasticityT300/%s/eps%s/itime%s\',\n"%(PHI[key],epsi,itime)
#			string[ inums ] = "\t\'2\':\'/../testRuns/glass%s\',\n"%(PHI[iphi])
	#
#				inums = lnums[ 2 ] - 1
#				string[ inums ] = "\targv2nd = \'itime=%s\\nindx=%s\'\n"%(itime*10000,key) #(PHI[iphi]*10000)
				string[ inums ] = "\targv2nd = \'itime=%s\\nitime0=%s\\nindx=%s\'\n"%(1000000,itime*nevery,6) #(PHI[iphi]*10000)

				sfile=open('junk%s.py'%count,'w');sfile.writelines(string);sfile.close()
				os.system( 'python3 junk%s.py'%count )
				os.system( 'rm junk%s.py'%count )
				count += 1
