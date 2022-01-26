import os

filee={1:'negative_mu',2:'mapCxyxy'}[1]
for i in [0,50,101,120,195]: #range(0,208,8):
	print('ElasticityT300/Co5Cr2Fe40Mn27Ni26/itime%i/Run0/%s.%i.png'%(i,filee,i*10000))
#	os.system('cp ElasticityT300/Co5Cr2Fe40Mn27Ni26/itime%i/Run0/%s.%s.png %s.%s.png'%(i,filee,i*10000,filee,i))
