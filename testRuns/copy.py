import os

dir='/mnt/home/share/nomaten/Kamran/glass'
n=3
write=['Co5Cr5Fe40Mn27Ni26', 'CoCrFeMn', 'CoNiCrFe', 'CoNiCrFeMn', 'CoNiFe', 'FeNi'][1]

for i in xrange(n):
	print 'irun=%s'%i
	os.system('mkdir -p %s/Run%s'%(write,i))
	os.system('cp %s/%s/calculated/%s_glass-300-100000000.dump %s/Run%s/%s_glass.dump'%(dir,write,i+1,write,i,write))
	os.system('cp %s/%s/calculated/%s_glass-300-100000000.txt  %s/Run%s/%s.txt'%(dir,write,i+1,write,i,write))
	os.system('cp %s/%s/calculated/%s_glass_300.data 		   %s/Run%s/%s_glass.data'%(dir,write,i+1,write,i,write))
