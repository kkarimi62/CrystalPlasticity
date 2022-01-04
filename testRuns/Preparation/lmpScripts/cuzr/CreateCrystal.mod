#--- comment this block if atoms are parsed from a data file
#--- Create box
variable 	a      			equal   ${cutoff}   #--- lattice constant
variable    volume          equal   ${natoms}*${a}^3.0/4.0 #--- natom * vol. of the voronoi cell
variable    lx              equal   floor(${volume}^(1.0/3.0)/${a})

#--- define crystal structure and lattice constant a0
#--- define direction vectors, i.e., set x=[100], y=[010], z=[001] and origin point.
lattice    fcc ${a} orient    x 1 0 0 orient y 0 1 0 orient z 0 0 1 &   
           origin 0.1 0.1 0.1
region    		mybox block 0 ${lx} 0 ${lx} 0 ${lx}   ## define box sizes along x, y, z (in the unit of a0)
create_box      2 mybox              ## create the simulation box, allowing a max of three species
create_atoms    1 box               ## create type-1 metal atoms in the box
#

#--- Need to set mass to something, just to satisfy LAMMPS
#mass 1 1.0e-20
mass            1  91.22 
mass		2  63.54 

#

#--- ternary alloys
#--- get number of atoms
variable	natom	equal	"atoms"
variable	natom2	equal	floor(${natom}*0.36)
#variable	natom3	equal	floor(${natom2})
# 
group kind1 type 1
set group kind1 type/subset 2 ${natom2} 12345
#
group kind2 type 2
group remain1 subtract all  kind2 
#
#set group remain1 type/subset 3 ${natom3} 74
#group kind3 type 3
#---- end of block

