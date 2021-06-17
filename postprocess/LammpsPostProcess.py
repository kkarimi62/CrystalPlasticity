#--- import system modules
import traceback
import sys
import numpy as np
import pandas as pd
import pdb

#--- utility funcs
def ConcatAttr( obj, attrs ):
    indx = 0
    existingAttrs = []
    for attr in attrs:
        if attr in dir(obj):
            existingAttrs += [ attr ]
            if indx == 0:
                XYZ_shifted = obj.__dict__[attr]
            else:
                XYZ_shifted = np.c_[XYZ_shifted,obj.__dict__[attr]]
            indx += 1
    return XYZ_shifted, existingAttrs

def shiftBeta( xyzDimensionless_j, diff ):
    indices = diff > 0.5
    beta0_j = xyzDimensionless_j - 1.0 * indices

    indices = diff < -0.5
    beta0_j += 1.0 * indices
    
    return beta0_j



############################################################
#######  class ReadDumpFile Reads LAMMPS dump files 
############################################################    
class ReadDumpFile:
    def __init__( self, path ):
        self.path = path
        self.coord_atoms_broken = {}
        self.BoxBounds = {}
    
    def GetCords( self, ncount = 1 ):
        slist = open( self.path )    
        count = 0
        try:
            while True and count <= ncount:
                sarr, cell_vector, itime, cols  = self.GetCordsTimeStep( slist ) #--- get coord

                #--- insert in a data frame
                self.coord_atoms_broken[ itime ] = pd.DataFrame( np.c_[sarr].astype('float'), columns = cols )

                #--- cast id and type to 'int'
                self.coord_atoms_broken[ itime ]['id'] = map(int,self.coord_atoms_broken[ itime ]['id'].tolist())[:]
                self.coord_atoms_broken[ itime ]['type'] = map(int,self.coord_atoms_broken[ itime ]['type'].tolist())[:]

                #--- sort
                self.coord_atoms_broken[ itime ].sort_values( by = 'id', inplace = True )

                #--- reset index
                self.coord_atoms_broken[ itime ].reset_index( drop=True, inplace=True )

                self.BoxBounds[ itime ] = cell_vector

                count += 1
        except:
#            traceback.print_exc()
            pass

    
    def GetCordsTimeStep(self, slist):
        slist.readline()
        itime = int( slist.readline().split()[0] )

        [slist.readline() for i in xrange(1)]
        nrows = int(slist.readline().split()[0])

        [slist.readline() for i in xrange(1)]

        CellVector = np.array([slist.readline().split() for i in xrange( 3 )])

        cols = slist.readline().split()[2:]

        return np.array([slist.readline().split() for i in xrange( nrows )]), CellVector, itime, cols

############################################################
#######  class with atom-related attributes 
############################################################    
class Atoms:
    def __init__(self,**kwargs):
#        print 'hello from Atoms const' 
        if 'x' in kwargs:
            self.x = kwargs['x']
        if 'y' in kwargs:
            self.y = kwargs['y']
        if 'z' in kwargs:
            self.z = kwargs['z']
        if 'xu' in kwargs:
            self.x = kwargs['xu']
        if 'yu' in kwargs:
            self.y = kwargs['yu']
        if 'zu' in kwargs:
            self.z = kwargs['zu']
        if 'id' in kwargs:
            self.id = kwargs['id']
        if 'type' in kwargs:
            self.type = kwargs['type']
        if 'xm' in kwargs:
            self.xm = kwargs['xm']
        if 'ym' in kwargs:
            self.ym = kwargs['ym']
        if 'zm' in kwargs:
            self.zm = kwargs['zm']
            
############################################################
#######  class with simulation cell attributes 
############################################################    
class Box:
    def __init__( self, **kwargs ):
		if 'BoxBounds' in kwargs:
			self.BoxBounds = kwargs['BoxBounds']
		if 'CellOrigin' in kwargs:
			self.CellOrigin = kwargs['CellOrigin']
		if 'CellVector' in kwargs:
			self.CellVector = kwargs['CellVector']
	#--- 2nd constructor
     
    def BasisVectors( self, **kwargs ):
		#     CellVector[0] = np.c_[CellVector[0],['0.0','0.0','0.0']] #--- ref. state
		if 'AddMissing' in kwargs:
			extraColumn = kwargs['AddMissing']
			if not self.BoxBounds.shape == (3,3):
				self.BoxBounds = np.c_[self.BoxBounds,extraColumn]
#				pdb.set_trace()
				print 'BoxBounds.shape=%s,%s is not (3,3)!'%(self.BoxBounds.shape)
				print 'add %s!'%(extraColumn)

		(xlo, xhi, xy) = map( float, self.BoxBounds[ 0 ] ) #--- xlo, xhi, xy
		lx = xhi - xlo - xy
		CellVector0 = np.array( [ lx, 0.0, 0.0 ] )

		(ylo, yhi, junk) =  map( float, self.BoxBounds[ 1 ] ) #--- ylo, yhi, xy
		ly = yhi - ylo
		a1 = np.array( [ 0.0, ly, 0.0 ] )
		CellVector1 = CellVector0 * ( xy / lx ) + a1

		(zlo, zhi, junk) =  map( float, self.BoxBounds[ 2 ] ) #--- zlo, zhi, xy
		lz = zhi - zlo
		CellVector2 = np.array( [ 0.0, 0.0, lz ] )

		self.CellOrigin = np.array( [ xlo, ylo, zlo ] )
		self.CellVector = np.c_[ CellVector0, CellVector1, CellVector2 ] 
        

class Wrap():
###########################################################
###### Wrap atom positions 
###########################################################    
    def __init__(self, atoms, box ):
        self.x = atoms.x
        self.y = atoms.y
        self.z = atoms.z
        self.CellVector = box.CellVector
        self.CellOrigin   = box.CellOrigin
  
    def GetDimensionlessCords( self ):
    ############################################################
    ####### compute dimensionless coords of atoms given 
    ####### corresponding cartesian coords
    ############################################################
        xyz_centered = np.c_[self.x,self.y,self.z] - self.CellOrigin
        self.beta = np.matmul( np.linalg.inv(self.CellVector), xyz_centered.T).T

    def WrapCoord( self ):
    ###########################################################
    ###### Wrap atom positions 
    ###########################################################    
        self.GetDimensionlessCords() #--- updates self.beta
        self.beta = self.beta % 1.0
        self.GetXYZ() #--- updates self.xyz
        indices = self.isInside()
        assert np.all( indices ), 'not all atoms are inside!'

    def GetXYZ( self ):
    ############################################################
    ####### compute coords of atoms given 
    ####### corresponding dimensionless coords
    ############################################################
        XYZ_centered = np.matmul( self.CellVector, self.beta.T ).T #--- xyz in reference state
        xyz = XYZ_centered + self.CellOrigin
        self.x = xyz[:,0]
        self.y = xyz[:,1]
        self.z = xyz[:,2]
        
    def isInside( self ):
        self.GetDimensionlessCords()
        #--- filter: only inside the expanded box 
        beta0 = self.beta[:,0]
        beta1 = self.beta[:,1]
        beta2 = self.beta[:,2]
        index0 = np.all([beta0 >= 0.0, beta0 < 1.0], axis=0 )
        index1 = np.all([beta1 >= 0.0, beta1 < 1.0], axis=0 )
        index2 = np.all([beta2 >= 0.0, beta2 < 1.0], axis=0 )
        return np.all([index0,index1,index2],axis=0)   

    def Set( self, atoms ):
        atoms.x = self.x
        atoms.y = self.y
        atoms.z = self.z

class Map(Wrap):
############################################################
####### map atoms within a tilted box to an orthogonal one
############################################################
    def __init__( self, atoms, box ):
        Wrap.__init__( self, atoms, box ) #--- call parent's constructor                  
    def ChangeBasis( self ):
        self.GetDimensionlessCords() #--- beta0, beta1, beta2
        #--- shift (beta0 + beta1 * dx / |b0| > 1) by - b0
        #--- [b0,b1]^{-1}*(x=|b0|,y)
        b2 = self.CellVector[:,2]
        b1 = self.CellVector[:,1]
        b0 = self.CellVector[:,0]
        norm_b0 = np.dot(b0,b0)**0.5
        shift0 = np.dot(b1, b0 / norm_b0 )
        #
        indices_shiftTrue = (self.beta[:,0]) + (self.beta[:,1])*shift0/norm_b0 >= 1.0
        n = len( indices_shiftTrue )
        shift_matrix =  np.array(indices_shiftTrue).reshape((n,1)) * b0
        #
        self.xm = np.c_[self.x,self.y,self.z] - shift_matrix
        
    def Set( self, atoms ):        
        sdict = {'xm':self.xm[:,0], 'ym':self.xm[:,1], 'zm':self.xm[:,2]}
        atoms.__init__(**sdict)
        

class Copy( Atoms, Wrap ):
############################################################
####### add replicals of the center simulation box
############################################################
    def __init__( self, atoms, box ):
        Atoms.__init__( self, **atoms.__dict__ ) #--- call parent's constructor
        Wrap.__init__( self, atoms, box )
        
    def FullCopies( self ):
        #--- concat. coordinates
        XYZ_shifted, attr0 = ConcatAttr( self, ['x','y','z','xu','yu','zu','xm','ym','zm'] )
        xyz_original = XYZ_shifted.copy()
        assert XYZ_shifted.shape[1] % 3 == 0, 'shifted coordinates must be integer multiple of 3!'
        #--- add other attributes to xyz 
#        attr_list = list(coord_atoms_broken.keys()) 

    #     #--- exclude xyz
        ID_TYPE_shifted, attr1 = ConcatAttr( self, ['id','type'])
        id_type_original = ID_TYPE_shifted.copy()
		
        #--- cell copies
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                for k in [-1,0,1]:
        #            print i,j,k
                    if i == j == k == 0:
                        continue
                    #
                    total_shift = np.matmul( self.CellVector, np.array([i,j,k]) ) #--- shape (3,)
                    if XYZ_shifted.shape[1] / 3 == 2:
                        total_shift = np.concatenate([total_shift,total_shift],axis=0) #--- shape: (6,)
                    if XYZ_shifted.shape[1] / 3 == 3:
                        total_shift = np.concatenate([total_shift,total_shift, total_shift],axis=0) #--- shape: (9,)
    #                assert total_shift.shape[ 0 ] == 9
                    #--- add shift
                    xyz_shifted = xyz_original + total_shift
                    #--- append
                    XYZ_shifted = np.concatenate( ( XYZ_shifted, xyz_shifted ), axis = 0 ) #--- bottleneck!!!!!
                    ID_TYPE_shifted = np.concatenate( ( ID_TYPE_shifted, id_type_original ), axis = 0 )
        #--- store full copies in a df
        assert ID_TYPE_shifted.shape[1] == len(attr1)
        assert XYZ_shifted.shape[1] == len(attr0)
        self.df = pd.DataFrame(np.c_[ID_TYPE_shifted,XYZ_shifted],columns=attr1 + attr0)

    def Get( self ):
        return Atoms( **self.df.to_dict(orient = 'list') ) #--- return atom object


    def Expand( self, ev = 0.1):
    ############################################################
    ####### Get atoms inside an expanded box
    ############################################################    
        assert 'xm' in dir(self) and 'ym' in dir(self) and 'zm' in dir(self), 'mapped coordinates are needed!'
        self.FullCopies() #--- full copies
        atomsCopied = self.Get() #--- new atom object

        indices =  self.isInsideExpanded(np.c_[atomsCopied.xm, atomsCopied.ym, atomsCopied.zm], #--- give mapped coordinates
                                         ev = ev,)
        #--- filter!!
        self.df = pd.DataFrame(atomsCopied.__dict__,)[indices]
        
    def isInsideExpanded( self, xyz, ev = 0.2 ):
    ############################################################
    ####### Get atoms inside an expanded box
    ############################################################    
        #--- tensor associated with dilation
        strainTensor = np.array([[ev/3.0,0,0],
                                [0,ev/3.0,0],
                                [0,0,ev/3.0]])
        #
        identityMat = np.array([[1,0,0],[0,1,0],[0,0,1]])
        #
        CellVectorOrtho, VectorNorm = GetOrthogonalBasis( self.CellVector ) #GetOrthogonalBasis?

        #--- extend diagonal
        rvect = -np.matmul( CellVectorOrtho, np.array([0.5,0.5,0.5]))
        CellOrigin_expanded = self.CellOrigin + np.matmul(strainTensor,rvect)

        #--- extend basis vectors
        CellVector_expanded = np.matmul( identityMat + strainTensor, CellVectorOrtho )

        
        #--- new atom, box, and wrap object
        atoms = Atoms(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2])
        box = Box( CellOrigin = CellOrigin_expanded, CellVector = CellVector_expanded )
        wrap = Wrap( atoms, box )
        
        return wrap.isInside()

class Compute( Atoms, Box ):
############################################################
####### base class for computing atom properties
############################################################
    def __init__( self, atoms, box ): #--- call parent constructor 
        Atoms.__init__(self, **atoms.__dict__) 
        Box.__init__(self, CellOrigin = box.CellOrigin, CellVector = box.CellVector )


class ComputeDisp( Compute, Wrap ):
############################################################
####### compute atoms' displacements given their
####### reference state
############################################################
    def __init__( self, atoms, box, atoms0, box0 ): #--- '0' denotes ref. state
        Compute.__init__( self, atoms, box ) #--- call parent's constructor        
        Wrap.__init__( self, atoms, box )
        self.atoms0 = atoms0
        self.box0 = box0
            
    def SetUnwrapped( self ):
        #--- only if unwrapped coords are available
        assert 'xu' in dir(self) and 'yu' in dir(self) and 'zu' in dir(self) and\
               'xu' in dir(self.atoms0) and 'yu' in dir(self.atoms0) and 'zu' in dir(self.atoms0),\
               'unwrapped coordinates are needed!'
       #--- displacement: r^{unwrpd}_j - r^{wrpd}_i
        disp = np.c_[self.xu,self.yu,self.zu] - np.c_[self.atoms0.xu,self.atoms0.yu,self.atoms0.zu]
        self.dx = disp[:,0]
        self.dy = disp[:,1]
        self.dz = disp[:,2]
        #--- 
        
    def SetWrapped( self ):
		self.EstimateUnwrappedCord()    
        #--- displacement: r^{unwrpd}_j - r^{wrpd}_i
		disp = np.c_[self.x,self.y,self.z] - np.c_[self.atoms0.x,self.atoms0.y,self.atoms0.z]
		self.dx = disp[:,0]
		self.dy = disp[:,1]
		self.dz = disp[:,2]
		print 'warning: attributes x, y, z are now unwrapped!'
        #--- 
 
    def EstimateUnwrappedCord( self ):
        #--- dimensionless cords
        self.GetDimensionlessCords() #--- current frame
        #
        wrap0 = Wrap(self.atoms0, self.box0)
        wrap0.GetDimensionlessCords() #--- reference frame
        
        #--- shift to get unwrapped cords
        diff = self.beta - wrap0.beta

        #--- new dimensionless cords
        beta0_j = shiftBeta( self.beta[:,0], diff[:,0])
        beta1_j = shiftBeta( self.beta[:,1], diff[:,1])
        beta2_j = shiftBeta( self.beta[:,2], diff[:,2])
        self.beta = np.c_[beta0_j,beta1_j,beta2_j]
        
        #--- unwrapped cords at deformed state
        self.GetXYZ()


if __name__ in '__main__':

	fileName = '/Users/Home/Desktop/Tmp/txt/git/CrystalPlasticity/BmgData/FeNi_glass.dump'
	myRDF = ReadDumpFile( fileName )
	myRDF.GetCords( ncount = sys.maxint )

	#
	myAtoms = Atoms( **myRDF.coord_atoms_broken[0].to_dict(orient='list') )
	#
	pdb.set_trace()
	myBox = Box( BoxBounds = myRDF.BoxBounds[0] )
	myBox.BasisVectors()


