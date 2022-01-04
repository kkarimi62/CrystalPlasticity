#--- compute stress
compute     peratom all stress/atom NULL
compute     pinit all reduce sum c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]
#
variable    press equal -(c_pinit[1]+c_pinit[2]+c_pinit[3])/(3*vol)
variable 	pxx0 atom -c_peratom[1]
variable 	pyy0 atom -c_peratom[2]
variable 	pzz0 atom -c_peratom[3]
variable 	pyz0 atom -c_peratom[6]
variable 	pxz0 atom -c_peratom[5]
variable 	pxy0 atom -c_peratom[4]

#--- store initial stress
#thermo_style	custom	step	v_press #v_sxx0[1] #	pxx #v_press

dump        1 all custom ${nevery} ${OUT_PATH}/dump_init.xyz id type x y z c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]
dump_modify 1 flush yes append yes format line "%d %d %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e"

#run	0

#undump		1

#variable	natom	equal	atoms
#variable iatom loop  ${natom} #--- six modes
#	label loop0_CompStrs0
#	#
#	variable	tmp		equal 	v_pxx0[${iatom}]
#	variable	sxx0_iatom${iatom}	equal ${tmp}
#	#
#	variable	tmp		equal 	v_pyy0[${iatom}]
#	variable	syy0_iatom${iatom}	equal ${tmp}
#	#
#	variable	tmp		equal 	v_pzz0[${iatom}]
#	variable	szz0_iatom${iatom}	equal ${tmp}
#	#
#	variable	tmp		equal 	v_pyz0[${iatom}]
#	variable	syz0_iatom${iatom}	equal ${tmp}
#	#
#	variable	tmp		equal 	v_pxz0[${iatom}]
#	variable	sxz0_iatom${iatom}	equal ${tmp}
#	#
#	variable	tmp		equal 	v_pxy0[${iatom}]
#	variable	sxy0_iatom${iatom}	equal ${tmp}
#	#
#
#next	iatom
#jump	${INC}/CompStrs0.mod loop0_CompStrs0


#uncompute	pinit
#uncompute	peratom


