import os
import numpy as np

filee={1:'negative_mu',2:'mapCxyxy'}[1]
glass={    '0':'FeNi',
           '1':'CoNiFe',
           '2':'CoNiCrFe',
           '3' :'CoCrFeMn',
           '4':'CoNiCrFeMn',
           '5':'Co5Cr2Fe40Mn27Ni26'
	  }['2']
irun=0
for i in [0,50,101,120,195]: #np.arange(0,200+1,1): #range(0,208,8): #[0,50,101,120,195]: #range(0,208,8):
	print('ElasticityT300/%s/itime%i/Run%s/%s.%i.png'%(glass,i,irun,filee,i*10000))
#	os.system('cp ElasticityT300/%s/itime%i/Run%s/%s.%s.png %s.%s.png'%(glass,i,irun,filee,i*10000,filee,i))
