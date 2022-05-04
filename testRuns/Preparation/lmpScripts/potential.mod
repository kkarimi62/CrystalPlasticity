# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potential
#pair_style	eam/alloy
#pair_coeff              * * ${PathEam}/ZrCu.lammps.eam Zr Cu
#
pair_style meam/c
pair_coeff * * ${PathEam}/library_CoNiCrFeMn.meam Co Ni Cr Fe Mn ${PathEam}/parameters_CoNiCrFeMn.meam Co Ni Cr Fe Mn

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		100
thermo_style custom step temp pe press xy vol #pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
