import pdb
import os

PHI=[1,2,4,6,8,10,12]
n=len(PHI)
for i in xrange(n):
	path = 'testCpuRunTime%s/Run0/log.lammps'%(i)
	os.system('grep %s %s>>junk.txt'%("\'Loop time of\'",path))


string = open('junk.txt').readlines()
try:
    for i in xrange(n):
        print PHI[i],string[i][13:string[i].find("on")-1]
except:
    pass
