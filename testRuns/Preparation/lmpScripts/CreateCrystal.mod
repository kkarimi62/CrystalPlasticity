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
create_box      ${ntype} mybox              ## create the simulation box, allowing a max of three species
create_atoms    1 box               ## create type-1 metal atoms in the box
#

#--- Need to set mass to something, just to satisfy LAMMPS
#mass    1  58.693    #--- Ni
#mass    2  58.933195 #--- Co
#mass    3  51.9961   #--- Cr

#--- unary alloys
#--- done

#--- binary alloys
#--- get number of atoms
variable	natom	equal	"atoms"
variable	natom2	equal	floor(${natom}/2)
#
group kind1 type 1
set group kind1 type/subset 2 ${natom2} ${seed0}
#---- end of block


#--- ternary alloys
#--- get number of atoms
#variable	natom	equal	"atoms"
#variable	natom2	equal	floor(${natom}/3)
#variable	natom3	equal	floor(${natom2})
#
#group kind1 type 1
#set group kind1 type/subset 2 ${natom2} ${seed0}
#
#group kind2 type 2
#group remain1 subtract all  kind2 
#
#set group remain1 type/subset 3 ${natom3} ${seed1}
#group kind3 type 3
#---- end of block

##--- 5-component alloy
##--- get number of atoms
#variable	natom	equal	"atoms"
#variable	natom2	equal	floor(${natom}/5)
#variable	natom3	equal	floor(${natom2})
#variable	natom4	equal	floor(${natom2})
#variable	natom5	equal	floor(${natom2})
##
#group kind1 type 1
#set group kind1 type/subset 2 ${natom2} ${seed0}
#group kind2 type 2
#
#group remain1 subtract all  kind2 
#set group remain1 type/subset 3 ${natom3} ${seed1}
#group kind3 type 3
#
#group remain2 subtract remain1  kind3 
#set group remain2 type/subset 4 ${natom4} ${seed2}
#group kind4 type 4
#
#group remain3 subtract remain2  kind4 
#set group remain3 type/subset 5 ${natom5} ${seed3}
#group kind5 type 5
##---- end of block

