import sys
import ovito
import ovito.modifiers as md
import numpy as np
import ovito.io as io #import import_file

from ovito.vis import Viewport, TachyonRenderer, RenderSettings

import math
import pdb

InputFile = sys.argv[1] 
OutputFile = sys.argv[2]
nevery = int(sys.argv[3])
AnalysisType = int(sys.argv[4]) #--- 0:CommonNeighborAnalysis 1:g(r) 2:d2min 3:voronoi analysis
if AnalysisType == 3:
    radii=list(map(float,sys.argv[5:]))

    
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
    d2min = md.AtomicStrainModifier(output_nonaffine_squared_displacements=True,
                                    eliminate_cell_deformation=True
                                   )
    d2min.reference.load(InputFile)
    pipeline.modifiers.append(d2min)

if AnalysisType == 3:
    # Set atomic radii (required for polydisperse Voronoi tessellation).
    #atypes = pipeline.source.particle_properties.particle_type.type_list
    type_property = pipeline.source.particle_properties.particle_type
#     print(radii)
    for t in type_property.type_list:
#         print(t.id)
        t.radius = radii[t.id-1]
    # Set up the Voronoi analysis modifier.
    voro = md.VoronoiAnalysisModifier(
                                    compute_indices = True,
                                    use_radii = True,
                                    edge_count = 9, # Length after which Voronoi index vectors are truncated
                                    edge_threshold = 0.1
                                    )
    pipeline.modifiers.append(voro)

    


for frame in range(0,pipeline.source.num_frames,nevery):
    # This loads the input data for the current frame and
    # evaluates the applied modifiers:
    pipeline.compute(frame)
    itime = pipeline.source.attributes['Timestep']
#    print(itime)
    
    if AnalysisType == 1:
        sfile.write('#ITIME\n%s\n'%itime)
        np.savetxt(sfile, cnm.rdf, header='r\tg(r)')

    # Access computed Voronoi indices as NumPy array.
    # This is an (N)x(edge_count) array.
#     if AnalysisType == 3:
#         voro_indices = pipeline.output.particle_properties['Voronoi Index'].array
    
#     pdb.set_trace()


if AnalysisType == 1:
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
    io.export_file( pipeline, OutputFile, "lammps_dump",\
                    columns = ["Particle Identifier", "Particle Type", "Position.X","Position.Y","Position.Z",\
                               "Nonaffine Squared Displacement"],
                     start_frame = 0,
                     end_frame = pipeline.source.num_frames,
                     every_nth_frame = nevery,
                     multiple_frames=True )

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
