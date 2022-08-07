import sys
import ovito
import ovito.modifiers as md
import numpy as np
import ovito.io as io #import import_file
from ovito.vis import Viewport, TachyonRenderer, RenderSettings
from ovito.data import CutoffNeighborFinder
import math
import pdb
import traceback


def GetNpairs(data, finder):        
    Npairs = 0
    for index in range(data.number_of_particles):
        for neigh in finder.find(index):
            Npairs += (index<neigh.index)
    return Npairs


def GetPairAttrs(data, neigh,iatom):
#    return list(map(lambda x:data.particle_properties.particle_identifier.array[x.index],neigh))
    return list(map(lambda x:(iatom,x.index,x.distance,x.delta[0],x.delta[1],x.delta[2], x.pbc_shift[0],x.pbc_shift[1],x.pbc_shift[2]),neigh))
        
        
#--- command-line args
InputFile = sys.argv[1] #--- input file
print('InputFile=',InputFile)
OutputFile = sys.argv[2] #--- output
print('OutputFile=',OutputFile)
nevery = int(sys.argv[3])
print('nevery',nevery)
AnalysisType = int(sys.argv[4]) #--- 0:CommonNeighborAnalysis 1:g(r) 2:d2min 3:voronoi analysis 4: neighborlist
if AnalysisType == 2 or AnalysisType == 5 or AnalysisType == 6: #--- d2min: reference file
    RefFile = sys.argv[5]

if AnalysisType == 3: #--- voronoi analysis
    radii=list(map(float,sys.argv[5:]))
if AnalysisType == 4: #--- neighbor list
    cutoff = float(sys.argv[5]) #--- cutoff distance

    
print('InputFile=',InputFile)
# Load input data and create a data pipeline.
pipeline = io.import_file('%s'%(InputFile), multiple_frames = True)
print('num_frames=',pipeline.source.num_frames)

# Calculate per-particle displacements with respect to initial simulation frame
if AnalysisType == 0:
    cna = md.CommonNeighborAnalysisModifier()
    pipeline.modifiers.append(cna)


#apply modifier
if AnalysisType == 1:
    cnm = md.CoordinationNumberModifier(cutoff = 10.0, number_of_bins = 200)
    pipeline.modifiers.append(cnm)
    sfile = open(OutputFile,'a')

if AnalysisType == 2:
    d2min = md.AtomicStrainModifier(
#                                    use_frame_offset = use_frame_offset, #True,
#                                    frame_offset = frame_offset, #-1,
#                                    reference_frame = reference_frame,
                                    output_nonaffine_squared_displacements=True,
                                    eliminate_cell_deformation=True,
                                   )
    d2min.reference.load(RefFile, multiple_frames = True)
    pipeline.modifiers.append(d2min)

if AnalysisType == 3:
    # Set atomic radii (required for polydisperse Voronoi tessellation).
    #atypes = pipeline.source.particle_properties.particle_type.type_list
    type_property = pipeline.source.particle_properties.particle_type
#     print(radii)
#    for t in type_property.type_list:
#         print(t.id)
#        t.radius = radii[t.id-1]
    # Set up the Voronoi analysis modifier.
    voro = md.VoronoiAnalysisModifier(
                                    compute_indices = True,
                                    use_radii = False, #True,
                                    edge_count = 9, # Length after which Voronoi index vectors are truncated
                                    edge_threshold = 0.1
                                    )
    pipeline.modifiers.append(voro)

#--- neighbor list
if AnalysisType == 4:
    sfile = open(OutputFile,'ab')

if AnalysisType == 5:
    strain = md.AtomicStrainModifier(
#                                    use_frame_offset = use_frame_offset, #True,
#                                    frame_offset = frame_offset, #-1,
#                                    reference_frame = reference_frame,
                                    output_strain_tensors=True,
                                    eliminate_cell_deformation=True,
                                   )
    strain.reference.load(RefFile, multiple_frames = True)
    pipeline.modifiers.append(strain)

if AnalysisType == 6:
    disp = md.CalculateDisplacementsModifier(
#                                    use_frame_offset = use_frame_offset, #True,
#                                    frame_offset = frame_offset, #-1,
#                                    reference_frame = reference_frame,
									affine_mapping = md.CalculateDisplacementsModifier.AffineMapping.ToReference
                                   )
    disp.reference.load(RefFile, multiple_frames = True)
    pipeline.modifiers.append(disp)
#if AnalysisType != 4:
#    frames = range(0,pipeline.source.num_frames,nevery)

#start_frame = 0
#for frame, indices in zip(frames,range(len(frames))):
for frame in range(0,pipeline.source.num_frames,nevery):
    # This loads the input data for the current frame and
    # evaluates the applied modifiers:
    print('frame=%s/%s'%(frame,pipeline.source.num_frames))
#    pipeline.compute(frame)
    try:
        data = pipeline.compute(frame)
    except:
        traceback.print_exc()
#        start_frame = frames[indices+1]
#        print('start_frame=',start_frame)
        continue
    itime = pipeline.source.attributes['Timestep']
#    print(itime)
    
    if AnalysisType == 1:
        sfile.write('#ITIME\n%s\n'%itime)
        np.savetxt(sfile, cnm.rdf, header='r\tg(r)')
        
        
    #--- compute neighbor list
    if AnalysisType == 4:
        type_property = pipeline.source.particle_properties.particle_type
        finder = CutoffNeighborFinder(cutoff, data)        
        neighList = list(map(lambda x: finder.find(x) , range(data.number_of_particles) ))
        zipp = zip(neighList,range(data.number_of_particles))
        pairij = np.concatenate(list(map(lambda x: GetPairAttrs( data, x[0],x[1] ), zipp))) #,dtype=object)
        #
        indexi = list(map(int,pairij[:,0]))
        indexj = list(map(int,pairij[:,1]))
        atomi_id=data.particle_properties.particle_identifier.array[indexi]
        atomj_id=data.particle_properties.particle_identifier.array[indexj]
        atomi_type = type_property.array[indexi]
        atomj_type = type_property.array[indexj]
        #
        sfile.write(b'ITIME: TIMESTEP\n%d\n'%itime)
        sfile.write(b'ITEM: NUMBER OF ATOMS\n%d\n'%(len(indexi)))
        sfile.write(b'ITEM: BOX BOUNDS xy xz yz pp pp pp\n0.0\t0.0\t0.0\n0.0\t0.0\t0.0\n0.0\t0.0\t0.0\n')
        sfile.write(b'ITEM: ATOMS id\ttype\tJ\tJtype\tDIST\tDX\tDY\tDZ\tPBC_SHIFT_X\tPBC_SHIFT_Y\tPBC_SHIFT_Z\n')
#        pdb.set_trace()
        np.savetxt(sfile,np.c_[ atomi_id, atomi_type, atomj_id, atomj_type, pairij[:,2:]],
                   fmt='%i %i %i %i %7.6e %7.6e %7.6e %7.6e %i %i %i' )
                    
#         for index in range(data.number_of_particles):
#             atomi_id = data.particle_properties.particle_identifier.array[index]
#             atomi_type = type_property.array[index]
#             if atomi_id == 1 or atomi_id == 2: print("Neighbors of particle %i:" % atomi_id)
# #            pdb.set_trace()
            # Iterate over the neighbors of the current particle:
#             for neigh in finder.find(index):
#                 atomj_id = data.particle_properties.particle_identifier.array[neigh.index]
#                 atomj_type = type_property.array[neigh.index]
#                 if atomi_id == 1 or atomi_id == 2: print("%i %i:" %(atomi_id,atomj_id))
#                 if atomi_id < atomj_id:
# #                print(neigh.index, neigh.distance, neigh.delta, neigh.pbc_shift)
#                     sfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(atomi_id,atomi_type,atomj_id, atomj_type, neigh.distance, neigh.delta[0],neigh.delta[1],neigh.delta[2], neigh.pbc_shift[0],neigh.pbc_shift[1],neigh.pbc_shift[2]))

    # Access computed Voronoi indices as NumPy array.
    # This is an (N)x(edge_count) array.
#     if AnalysisType == 3:
#         voro_indices = pipeline.output.particle_properties['Voronoi Index'].array
    
#    pdb.set_trace()


if AnalysisType == 1 or AnalysisType == 4:
    sfile.close()
    
#--- export data
if AnalysisType == 0:
    io.export_file( pipeline, OutputFile, "lammps_dump",\
                    columns = ["Particle Identifier", "Particle Type", "Position.X","Position.Y","Position.Z",\
                               "Structure Type"],
                     start_frame = 0,
                     end_frame = pipeline.source.num_frames,
                     every_nth_frame = nevery,
                     multiple_frames=True )
if AnalysisType == 2:
#    print('start_frame=',start_frame)
#    print(help(io.export_file))
    io.export_file( pipeline, OutputFile, "lammps_dump",\
                    columns = ["Particle Identifier", "Particle Type", "Position.X","Position.Y","Position.Z",\
                               "Nonaffine Squared Displacement"],
                     start_frame = 0,
                     end_frame = pipeline.source.num_frames,
                     every_nth_frame = nevery,
                     multiple_frames=True 
                  )

if AnalysisType == 3: 
    io.export_file( pipeline, OutputFile, "lammps_dump",\
                    columns = ["Particle Identifier", "Particle Type", "Position.X","Position.Y","Position.Z",\
                               "Voronoi Index.0","Voronoi Index.1","Voronoi Index.2",\
                               "Voronoi Index.3","Voronoi Index.4","Voronoi Index.5",\
                               "Voronoi Index.6","Voronoi Index.7","Voronoi Index.8", "Atomic Volume"],
                     start_frame = 0,
                     end_frame = pipeline.source.num_frames,
                     every_nth_frame = nevery,
                    
                    multiple_frames=True 
                  )   

if AnalysisType == 5:
#    print('start_frame=',start_frame)
#    print(help(io.export_file))
    io.export_file( pipeline, OutputFile, "lammps_dump",\
                    columns = ["Particle Identifier", "Particle Type", "Position.X","Position.Y","Position.Z",\
                               "Strain Tensor.XX", "Strain Tensor.YY","Strain Tensor.ZZ", "Strain Tensor.XY",
							   "Strain Tensor.XZ", "Strain Tensor.YZ"
							  ],
                     start_frame = 0,
                     end_frame = pipeline.source.num_frames,
                     every_nth_frame = nevery,
                     multiple_frames=True 
                  )

if AnalysisType == 6:
#    print('start_frame=',start_frame)
#    print(help(io.export_file))
    io.export_file( pipeline, OutputFile, "lammps_dump",\
                    columns = ["Particle Identifier", "Particle Type", "Position.X","Position.Y","Position.Z",\
                               "Displacement.X", "Displacement.Y","Displacement.Z"],
                     start_frame = 0,
                     end_frame = pipeline.source.num_frames,
                     every_nth_frame = nevery,
                     multiple_frames=True 
                  )
# Export the computed RDF data to a text file.

'''
pipeline.dataset.anim.frames_per_second = 60
pipeline.add_to_scene()
vp = Viewport()

vp.type = Viewport.Type.PERSPECTIVE

#vp.camera_pos = (735.866,-725.04,1001.35)
vp.camera_pos = (118.188,-157.588,131.323)

#vp.camera_dir = (-0.49923, 0.66564, -0.5547)
vp.camera_dir = (-0.49923,0.66564,-0.5547) 

vp.fov = math.radians(35.0)

tachyon = TachyonRenderer() #shadows=False, direct_light_intensity=1.1)

rs = RenderSettings(size=(600,600), filename="image.mov",
#                   custom_range=(0,100),
                    everyNthFrame=1,
                    range = RenderSettings.Range.ANIMATION, #CUSTOM_INTERVAL, #RenderSettings.Range.ANIMATION,  
                    renderer=tachyon,
                    )

vp.render(rs)
'''
