import numpy as np
# 
def Generate( natom, ntypes, 
#			  (xlo,xhi), (ylo, yhi), (zlo, zhi),
			  rho,
			  title = 'data.txt',
			  **kwargs ): 

	L = (natom/rho) ** (1.0/3.0)
	xlo = ylo = zlo = 0.0
	xhi = yhi = zhi = L

	x, y, z  = np.random.rand(	natom, 3 ).T

	x *= ( xhi - xlo )
	x += xlo
	y *= ( yhi - ylo )
	y += ylo
	z *= ( zhi - zlo )
	z += zlo
	
	#--- assign types
	assert len(kwargs) == ntypes, 'provide number ration for each species!'
	types = np.zeros((natom,), dtype=int)
	atomi = 0
	for key in kwargs:
		ratio = kwargs[ key ]
		itype = key[-1]

		atomj = atomi + int( ratio * natom )
		types[atomi:atomj] = itype

		atomi = atomj
	assert atomj <= natom
	types[atomj:natom] = itype
	#--- write output
	sfile = open( title, 'w' )
	print >> sfile, '# LAMMPS data file'
	print >> sfile, '%s atoms'%natom
	print >> sfile, '%s atom types'%ntypes
	print >> sfile, '%s %s xlo xhi'%(xlo,xhi)
	print >> sfile, '%s %s ylo yhi'%(ylo,yhi)
	print >> sfile, '%s %s zlo zhi'%(zlo,zhi)
	print >> sfile
	print >> sfile, 'Masses'
	print >> sfile
	for i in xrange( ntypes ):
		print >> sfile, '%s %s'%(i+1,1.0)
	print >> sfile
	print >> sfile, 'Atoms  # atomic: id tyoe x y z'
	print >> sfile
	for iatom in xrange( natom ): 
		print >> sfile, '%s %s %s %s %s'%( iatom+1, types[ iatom ], x[ iatom ], y[ iatom ], z[ iatom ])
	sfile.close()
	
if __name__ == '__main__':

	Generate( 1000, 5, 
			 (0,1.0), (0, 1.0), (0.0, 1.0),
			 title = 'data.txt',
			 ratio1 = 0.05, ratio2 = 0.26, ratio3 = 0.02, ratio4 = 0.4, ratio5 = 0.27 )	
	
