import sys
import ovito
import ovito.modifiers as md
import numpy as np
import ovito.io as io #import import_file

from ovito.vis import Viewport, TachyonRenderer, RenderSettings
import os
import math
import pdb

path='.' #'/Users/Home/Desktop/Tmp/txt/scratch/10561052'
InputFile = sys.argv[1] 
OutputFile = sys.argv[2]
AnalysisType = int(sys.argv[3]) #--- 0:CommonNeighborAnalysis 1:CoordinationNumber

# Load input data and create a data pipeline.
pipeline = io.import_file('%s/%s'%(path,InputFile), multiple_frames = True)
print('num_frames=',pipeline.source.num_frames)

# Calculate per-particle displacements with respect to initial simulation frame
if AnalysisType == 0:
	cna = md.CommonNeighborAnalysisModifier()
	pipeline.modifiers.append(cna)

#apply modifier
if AnalysisType == 1:
	cnm = md.CoordinationNumberModifier(cutoff = 10.0, number_of_bins = 200)
	pipeline.modifiers.append(cnm)
	try:
		os.system('rm %s'%OutputFile)
	except:
		pass
	sfile = open(OutputFile,'a')
	
for frame in range(pipeline.source.num_frames):
#	pdb.set_trace()
    # This loads the input data for the current frame and
    # evaluates the applied modifiers:
	pipeline.compute(frame)
	itime = pipeline.source.attributes['Timestep']
	if AnalysisType == 1:
		sfile.write('#ITIME\n%s\n#r\tg(r)\n'%itime)
		for item in cnm.rdf:
			sfile.write('%s\t%s\n'%(item[0],item[1]))
#		np.savetxt(sfile, cnm.rdf, header="r g(r)")

if AnalysisType == 1:
	sfile.close()	
#--- export data
if AnalysisType == 0:
	io.export_file( pipeline, OutputFile, "lammps_dump",\
					columns = ["Particle Identifier", "Particle Type", "Structure Type"],
					multiple_frames=True )

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
