import pdb
import os

PHI=range(2,32+2,2)
n=len(PHI)
for i in xrange(n):
	path = 'testCpuRunTime%s/Run0/log.lammps'%(i)
	os.system('grep %s %s>>junk.txt'%("\'Loop time of\'",path))


string = open('junk.txt').readlines()
try:
    for i in xrange(n):
        print string[i][string[i].find("on")+3:string[i].find("procs")-1] ,string[i][13:string[i].find("on")-1]
except:
    pass
