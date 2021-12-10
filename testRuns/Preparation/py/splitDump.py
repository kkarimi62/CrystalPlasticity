def GetAtoms( filee, nevery = 1 ):
    lmpData = lp.ReadDumpFile( filee )
    lmpData.GetCords( ncount = sys.maxsize, 
                     columns = {'c_peratom[1]':'sxx','c_peratom[2]':'syy','c_peratom[3]':'szz',
                                                       'c_peratom[4]':'sxy','c_peratom[5]':'sxz','c_peratom[6]':'syz'}
                    )
    nn=len(lmpData.coord_atoms_broken.keys())
    itimee=list(lmpData.coord_atoms_broken.keys())[0:nn:nevery]

    #--- atom obj
    box0 = map(lambda x:lp.Box( BoxBounds = lmpData.BoxBounds[x], AddMissing = np.array([0.0,0.0,0.0] ) ), itimee ) #--- reference state
    atoms0 = map(lambda x: lp.Atoms( **lmpData.coord_atoms_broken[x].to_dict(orient='series')), itimee )

    return dict(zip(itimee,list(atoms0))), dict(zip(itimee,list(box0)))


#--- user modules
import sys
sys.path += ['/home/kamran.karimi1/Project/git/CrystalPlasticity/postprocess']
import LammpsPostProcess2nd as lp
import numpy as np
#import imp
#imp.reload(lp)

input_file = sys.argv[1]
mass = {'1':58.933,'2':58.690, '3':51.960, '4':55.847,'5':54.940}

atomi, boxi = GetAtoms( input_file, nevery = 1 ) #--- change it to 1
#---
count = 0
for key in atomi:
    wd = lp.WriteDataFile( atomi[key], boxi[key], mass )
    wd.Write('data.%s.txt'%count)
    count += 1
