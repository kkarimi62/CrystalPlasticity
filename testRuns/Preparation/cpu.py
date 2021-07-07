import pdb
import os

PHI=[1,2,4,6,8,10,12,14,16,18,20,22,24,26,28] #range(0,2000000+10000,4*10000) 
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
