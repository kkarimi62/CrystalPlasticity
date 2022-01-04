#--- comment if the box and atoms are created below!
if "${ParseData} == 1" then &
	"read_data	${DataFile}" & 
else &
	"include ${INC}/CreateCrystal.mod"

#---
#change_box	all	triclinic

