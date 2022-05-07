#--- compute stress
compute     peratom all stress/atom NULL
compute     p all reduce sum c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]
#
variable    press equal -(c_p[1]+c_p[2]+c_p[3])/(3*vol)
variable 	pxx1_${dir} atom -c_peratom[1]
variable 	pyy1_${dir} atom -c_peratom[2]
variable 	pzz1_${dir} atom -c_peratom[3]
variable 	pyz1_${dir} atom -c_peratom[6]
variable 	pxz1_${dir} atom -c_peratom[5]
variable 	pxy1_${dir} atom -c_peratom[4]

dump        2 all custom ${nevery} ${OUT_PATH}/${DumpFile}${dir}.xyz id type x y z c_peratom[1] c_peratom[2] c_peratom[3] c_peratom[4] c_peratom[5] c_peratom[6]
dump_modify 2 flush yes append yes format line "%d %d %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e"

#--- store initial stress
#thermo_style	custom	step	v_press
#run	0

#undump		1


#variable	pxx1_${dir}_atom0	equal v_sxx1_${dir}[1]
#variable	pxx1_${dir}_atom0	equal v_sxx1_${dir}[1]
#variable	pyy1_${dir}	atom ${pyy0}
#variable	pzz1_${dir}	atom ${pzz0}
#variable	pyz1_${dir}	atom ${pyz0}
#variable	pxz1_${dir}	atom ${pxz0}
#variable	pxy1_${dir}	atom ${pxy0}

#uncompute	p
#uncompute	peratom

