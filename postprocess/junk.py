import os
filee='cr_mu.png' #'map_d2min.png' #'pdf.pmg'
for i in [0,50,100]:
	os.system('cp FeNiT300Elasticity%s/Run1/%s figure.%s.png'%(i,filee,i))
