# NOTE: This script should not need to be
# modified. See in.elastic for more info.
#
# Find which reference length to use

if "${dir} == 1" then &
   "variable len0 equal ${lx0}" 
if "${dir} == 2" then &
   "variable len0 equal ${ly0}" 
if "${dir} == 3" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 4" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 5" then &
   "variable len0 equal ${lz0}" 
if "${dir} == 6" then &
   "variable len0 equal ${ly0}" 

# Reset box and simulation parameters

clear
box tilt large
read_restart restart.equil
include ${INC}/potential.mod
change_box all triclinic

#dump        1 all custom 1 dump.xyz id type x y z
#dump_modify 1 flush yes append yes
#run	0
#undump	1

# Negative deformation

variable delta equal -${up}*${len0}
variable deltaxy equal -${up}*xy
variable deltaxz equal -${up}*xz
variable deltayz equal -${up}*yz
if "${dir} == 1" then &
   "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"
if "${dir} == 2" then &
   "change_box all y delta 0 ${delta} yz delta ${deltayz} remap units box"
if "${dir} == 3" then &
   "change_box all z delta 0 ${delta} remap units box"
if "${dir} == 4" then &
   "change_box all yz delta ${delta} remap units box"
if "${dir} == 5" then &
   "change_box all xz delta ${delta} remap units box"
if "${dir} == 6" then &
   "change_box all xy delta ${delta} remap units box"



# Relax atoms positions

#velocity	all       create     $T     97364  mom    yes  rot  yes  dist  gaussian  # Dynamics starts
#fix			FixTherm        all        nvt temp   $T  $T  ${damp_t}  
#fix fixLangevin all langevin ${temp} ${temp} ${tdamp} ${seed}
#run ${nrun}

#minimize ${etol} ${ftol} ${maxiter} ${maxeval}

#dump        1 all custom 1 dump.xyz id type x y z
#dump_modify 1 flush yes append yes
#run	0

