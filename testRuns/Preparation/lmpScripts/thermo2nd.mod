#--- output thermodynamic variables
variable varStep equal step
#variable varPress equal press
variable varTemp equal temp
variable varTime equal v_varStep*${dt}
#variable varXy	 equal	xy
#variable varLy	 equal	ly
variable varVol	 equal	vol
#variable varPe	 equal	pe
#variable varPxx	 equal	pxx
#variable varPyy	 equal	pyy
#variable varPzz	 equal	pzz
#variable varPxy	 equal	pxy
#variable varSxy	 equal	-v_varPxy*${cfac}
#variable varSzz	 equal	-v_varPzz*${cfac}
#variable varExy	 equal	v_varXy/v_varLy		
#variable varEzz	 equal	v_varTime*${GammaDot}		
#variable ntherm  equal	ceil(${Nstep}/${nthermo})
#variable varn	 equal	v_ntherm
variable       swap_attempt equal f_4[1]
variable       swap_accept  equal f_4[2]
variable 	   var_swap_attempt equal v_swap_attempt
variable 	   var_swap_accept  equal v_swap_accept
#compute thermo_temp0 bulk temp
#variable varTemp equal c_thermo_temp0


fix extra all print 10 "${varStep} ${varTime} ${varTemp} ${varVol} ${var_swap_accept} ${var_swap_attempt}" screen no title "step time temp vol swap_accept swap_attempt" file thermo2nd.txt

thermo 10 #0
thermo_style custom step varTemp pe press vol v_swap_accept v_swap_attempt
thermo_modify norm no

