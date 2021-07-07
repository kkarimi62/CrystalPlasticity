
import os

PHI=[1,2,4,6,8,10,12]
n=len(PHI)
for i in xrange(n):
	path = 'testCpuRunTime%s/Run0/log.lammps'%(i)
	os.system('grep %s %s'%('Loop time of',path))

