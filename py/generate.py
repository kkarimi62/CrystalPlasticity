import numpy as np
import pandas as pd
import sys
import pdb
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
	print >> sfile, 'Atoms  # atomic: id type x y z'
	print >> sfile
	for iatom in xrange( natom ): 
		print >> sfile, '%s %s %s %s %s'%( iatom+1, types[ iatom ], x[ iatom ], y[ iatom ], z[ iatom ])
	sfile.close()

def GenerateDataFromDump( pathDump, pathData = '', 
			  title = 'data.txt',
			  **kwargs ): 
	coord, cell_vector = GetCords( pathDump )
        coord[0]=coord[coord.keys()[0]]
        cell_vector[0]=cell_vector[coord.keys()[0]]
#        pdb.set_trace()
	CellOrigin, CellVector = GetCellVector( cell_vector )
	
	(xlo, ylo, zlo) = CellOrigin[ 0 ]
	xhi = (CellOrigin[ 0 ] + CellVector[0][:,0])[ 0 ]
	yhi = (CellOrigin[ 0 ] + CellVector[0][:,1])[ 1 ]
	zhi = (CellOrigin[ 0 ] + CellVector[0][:,2])[ 2 ]

	x, y, z  = np.c_[coord[0][['x','y','z']]].T
        if 'vx' in coord[0].keys():
    	    vx, vy, vz  = np.c_[coord[0][['vx','vy','vz']]].T
        else:
            vx = vy = vz = np.zeros(len( coord[0] ))

	types = coord[0]['type'].to_list()
	ntypes = kwargs['ntype'] if 'ntype' in kwargs else len(set(types))
	
	natom = len( coord[0] )
	#--- write output
	sfile = open( title, 'w' )
	print >> sfile, '# LAMMPS data file'
	print >> sfile, '%s atoms'%natom
	print >> sfile, '%s atom types'%ntypes
	print >> sfile, '%s %s xlo xhi'%(xlo,xhi)
	print >> sfile, '%s %s ylo yhi'%(ylo,yhi)
	print >> sfile, '%s %s zlo zhi'%(zlo,zhi)
#	print >> sfile
#	print >> sfile, 'Masses'
#	print >> sfile
#	for i in xrange( ntypes ):
#		print >> sfile, '%s %s'%(i+1,1.0)
	print >> sfile
	print >> sfile, 'Atoms  # atomic: id type x y z'
	print >> sfile
	for iatom in xrange( natom ): 
		print >> sfile, '%s %s %s %s %s'%( iatom+1, types[ iatom ], x[ iatom ], y[ iatom ], z[ iatom ])

	print >> sfile
	print >> sfile, 'Velocities  # atomic: id  vx vy vz'
	print >> sfile
	for iatom in xrange( natom ): 
		print >> sfile, '%s %s %s %s'%( iatom+1, vx[ iatom ], vy[ iatom ], vz[ iatom ])
	sfile.close()
	
def GetCords( file_name,
              ncount = sys.maxint):
    slist = open( file_name )    
    coord_atoms_broken = {}
    CellVector = {}
    count = 0
    try:
        while True and count <= ncount:
            sarr, cell_vector, itime, cols  = GetCordsTimeStep( slist ) #--- get coord
            #--- insert in a data frame
            coord_atoms_broken[ itime ] = pd.DataFrame( sarr, columns = cols )
            
            #--- cast id to 'int'
            coord_atoms_broken[ itime ]['id'] = map(int,coord_atoms_broken[ itime ]['id'].tolist())[:]

            #--- sort
            coord_atoms_broken[ itime ].sort_values( by = 'id', inplace = True )
            
            #--- reset index
            coord_atoms_broken[ itime ].reset_index( drop=True, inplace=True )
            
            CellVector[ itime ] = cell_vector
            
            count += 1
    except:
#        traceback.print_exc()
        pass

    
    return coord_atoms_broken, CellVector

def GetOrthogonalBasis( CellVector ): 
    a0 = CellVector[:,0]
    l0 = np.dot(a0,a0)**0.5
    n0 = a0 / l0 
    #
    a1 = CellVector[:,1]
    a1_perp = a1 - np.dot( a1, n0 ) * n0
    l1 = np.dot( a1_perp, a1_perp) ** 0.5
    #
    a2 = CellVector[:,2]
    l2 = np.dot(a2,a2)**0.5
    
    return np.c_[a0, a1_perp, a2], [l0, l1, l2]

def GetCordsTimeStep(slist):
    slist.readline()
    itime = int( slist.readline().split()[0] )
    
    [slist.readline() for i in xrange(1)]
    nrows = int(slist.readline().split()[0])

    [slist.readline() for i in xrange(1)]
    
    CellVector = np.array([slist.readline().split() for i in xrange( 3 )])
    
    cols = slist.readline().split()[2:]
    
    return np.array([slist.readline().split() for i in xrange( nrows )]), CellVector, itime, cols

def GetCellVector( CellVector ):
#    pdb.set_trace()
    CellVector[0] = np.c_[CellVector[0],['0.0','0.0','0.0']] #--- ref. state

    CellOrigin = {}
    cell_vector = {}
    
    for itime in CellVector:
        row0 = map( float, CellVector[ itime ][ 0 ] )
        l0 = row0[ 1 ] - row0[ 0 ]# - row0[ 2 ]
        CellVector0 = np.array( [ l0, 0.0, 0.0 ] )
        dx = 0.0 #row0[ 2 ]
        
        row1 =  map( float, CellVector[ itime ][ 1 ] )
        l1 = row1[ 1 ] - row1[ 0 ]
        a1 = np.array( [ 0.0, l1, 0.0 ] )
        CellVector1 = CellVector0 * ( dx / l0 ) + a1
        
        row2 =  map( float, CellVector[ itime ][ 2 ] )
        l2 = row2[ 1 ] - row2[ 0 ]
        CellVector2 = np.array( [ 0.0, 0.0, l2 ] )
        
        CellOrigin[ itime ] = np.array( [ row0[ 0 ], row1[ 0 ], row2[ 0 ] ] )
        cell_vector[ itime ] = np.c_[ CellVector0, CellVector1, CellVector2 ] 

    return CellOrigin, cell_vector


if __name__ == '__main__':

#	Generate( 1000, 5, 
#			 (0,1.0), (0, 1.0), (0.0, 1.0),
#			 title = 'data.txt',
#			 ratio1 = 0.05, ratio2 = 0.26, ratio3 = 0.02, ratio4 = 0.4, ratio5 = 0.27 )	

        pathh = '/home/kamran.karimi1/Project/git/CrystalPlasticity/testRuns/Preparation/test11thMelt2ndNatom50kParallel/Run0/AsQuenched.dump'
	GenerateDataFromDump( pathh, pathh, ntype = 5, 
			  title = 'junk.txt') 
	
