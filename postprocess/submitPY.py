if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    #---
    lnums = [ 40, 49 ]
    string=open('postprocess_ncbj_slurm.py').readlines() #--- python script
    #---
    PHI ={ 
#            '0':'FeNi',
#            '1':'CoNiFe',
#           '2':'CoNiCrFe',
#           '3' :'CoCrFeMn',
#            '4':'CoNiCrFeMn',
            '5':'Co5Cr2Fe40Mn27Ni26'
#			'6':'cuzr',
        }

    temps={
#        0:298,
        1:423,
        2:523,
        3:623,
        4:723
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
    times=np.arange(0,200+1,2)  #--- run GetFrames(lmpData,times=lmpData.coord_atoms_broken.keys())
    nevery = 10000 #--- run list(lmpData.coord_atoms_broken.keys())[1]

    count = 0
    for key_t in temps:
        temp = temps[key_t]
# 		for epsi in EPS:
# 			for itime in times:
            #---	
                inums = lnums[ 0 ] - 1
#				string[ inums ] = "\t\'6\':\'MlTrain/%s/itime%s',\n" % (PHI[key],itime) #--- change job name
                string[ inums ] = "\t\'9\':\'amorphousAlumina/temp%s\',\n" % (key_t) #--- change job name
#			string[ inums ] = "\t\'3\':\'PairCrltnT300/%s',\n" % (PHI[iphi]) #--- change job name
    #		string[ inums ] = "\t\'4\':\'VorAnlT300/%s',\n" % (PHI[iphi]) #--- change job name
    #		string[ inums ] = "\t\'5\':\'D2minAnalysisT300/%s',\n" % (PHI[iphi]) #--- change job name
        #---	densities
                inums = lnums[ 1 ] - 1
                string[ inums ] = "\t\'9\':\'/../testRuns/glassAlumina/temp%s\',\n"%(key_t)
#			string[ inums ] = "\t\'2\':\'/../testRuns/glass%s\',\n"%(PHI[iphi])
    #
#				inums = lnums[ 2 ] - 1
#				string[ inums ] = "\targv2nd = \'itime=%s\\nindx=%s\'\n"%(itime*10000,key) #(PHI[iphi]*10000)
#				string[ inums ] = "\targv2nd = \'itime=%s\\nitime0=%s\\nindx=%s\'\n"%(1000000,itime*nevery,6) #(PHI[iphi]*10000)
# 				inums = lnums[ 2 ] - 1
# 				string[ inums ] = "\tconfParser.set(\'parameters\',\'itime0\',\'%s\')\n"%(0) 
# 				inums = lnums[ 3 ] - 1
# 				string[ inums ] = "\tconfParser.set(\'parameters\',\'itime\',\'%s\')\n"%(itime*nevery) 
# 				inums = lnums[ 4 ] - 1
# 				string[ inums ] = "\tconfParser.set(\'input files\',\'fileIndex\',\'%s\')\n"%(key) 

                sfile=open('junk%s.py'%count,'w');sfile.writelines(string);sfile.close()
                os.system( 'python3 junk%s.py'%count )
                os.system( 'rm junk%s.py'%count )
                count += 1
