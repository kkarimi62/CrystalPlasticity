if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    #---
    lnums = [ 40, 49 ]
    string=open('postprocess_ncbj_slurm.py').readlines() #--- python script
    #---

    temps={
#        0:298,
        1:423,
        2:523,
        3:623,
        4:723
    }




    count = 0
    for key_t in temps:
        temp = temps[key_t]
                inums = lnums[ 0 ] - 1
                string[ inums ] = "\t\'9\':\'amorphousAlumina/temp%s\',\n" % (key_t) #--- change job name

                inums = lnums[ 1 ] - 1
                string[ inums ] = "\t\'9\':\'/../testRuns/glassAlumina/temp%s\',\n"%(key_t)

                sfile=open('junk%s.py'%count,'w');sfile.writelines(string);sfile.close()
                os.system( 'python3 junk%s.py'%count )
                os.system( 'rm junk%s.py'%count )
                count += 1
