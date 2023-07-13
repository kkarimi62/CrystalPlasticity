def makeOAR( EXEC_DIR, node, core, time ):
    someFile = open( 'oarScript.sh', 'w' )
    print >> someFile, '#!/bin/bash\n'
    print >> someFile, 'EXEC_DIR=%s\n' %( EXEC_DIR )
    print >> someFile, 'MEAM_library_DIR=%s\n' %( MEAM_library_DIR )
#	print >> someFile, 'module load mpich/3.2.1-gnu\n'
    print >> someFile, 'module load openmpi/4.0.2-gnu730\n'

    #--- run python script 
#	 print >> someFile, "$EXEC_DIR/%s < in.txt -var OUT_PATH %s -var MEAM_library_DIR %s"%( EXEC, OUT_PATH, MEAM_library_DIR )
#	cutoff = 1.0 / rho ** (1.0/3.0)
    for script,var,indx, execc in zip(Pipeline,Variables,range(100),EXEC):
        if execc == 'lmp': #_mpi' or EXEC == 'lmp_serial':
            print >> someFile, "mpirun --oversubscribe -np %s $EXEC_DIR/lmp_mpi < %s -echo screen -var OUT_PATH %s -var PathEam %s -var INC %s %s"%(nThreads*nNode, script, OUT_PATH, MEAM_library_DIR, SCRPT_DIR, var)
        elif execc == 'py':
            print >> someFile, "python3 %s %s"%(script, var)

    someFile.close()										  


if __name__ == '__main__':
    import os
    import numpy as np

    runs	 = [0] #,1,2]
    #
    nThreads = 2
    nNode	 = 1
    #
    jobname  = {
                1:'CuZrNatom32KT300Tdot1E-3Sheared',
                2:'CuZrNatom32KT300Tdot1E-1Elasticity',
                4:'ElasticityT300/Co5Cr2Fe40Mn27Ni26/itime0',
                5:'annealing',
               }[5]
    sourcePath = os.getcwd() +\
                {	
                    1:'/../postprocess/NiCoCrNatom1K',
                    2:'/CuZrNatom32KT300Tdot1E-1Sheared',
                    3:'/../glassCo5Cr2Fe40Mn27Ni26',
                    4:'/junk',
                }[4] #--- must be different than sourcePath
        #
    sourceFiles = { 0:False,
                    1:['Equilibrated_300.dat'],
                    2:['data.txt','ScriptGroup.txt'],
                    3:['data.0.txt','Co5Cr2Fe40Mn27Ni26_glass.dump','Co5Cr2Fe40Mn27Ni26.txt'], 
                    4:['data_minimized.txt'],
                    5:['data_init.txt','ScriptGroup.0.txt'], #--- only one partition! for multiple ones, use 'submit.py'
                    6:['data.0.txt','dumpSheared.xyz'], 
                 }[0] #--- to be copied from the above directory
    #
    EXEC_DIR = '/home/kamran.karimi1/Project/git/lammps2nd/lammps/src' #--- path for executable file
    #
    MEAM_library_DIR='/home/kamran.karimi1/Project/git/lammps2nd/lammps/potentials'
    #
    SCRPT_DIR = os.getcwd()+'/lmpScripts' #/'+{1:'cuzr'}[1]
    #
    SCRATCH = True
    OUT_PATH = '.'
    if SCRATCH:
        OUT_PATH = '/scratch/${SLURM_JOB_ID}'
    #--- py script must have a key of type str!
    LmpScript = {	0:'in.melt',
                    1:'relax.in', 
                    2:'relaxWalls.in', 
                    7:'in.Thermalization', 
                    4:'in.vsgc', 
                    5:'in.minimization', 
                    6:'in.shearDispTemp', 
                    8:'in.shearLoadTemp',
                    9:'in.elastic',
                    10:'in.elasticTemp',
                    101:'in.elasticTemp',
                    11:'in.oscillatoryShear',
                    13:'in.swap_GB4',
                    'p0':'partition.py',
                    'p1':'WriteDump.py',
                    'p2':'DislocateEdge.py',
                } 
    #
    Variable = {
                0:'  -var ParseData 1 -var DataFile equilibrated.txt -var tstart 300.0 -var tstop 2000.0 -var TdotMelt 100.0 -var TdotQuench 0.001 -var Pinit 1.0132 -var nevery 10000  -var DumpFile dumpInit.xyz -var WriteData data_quenched.dat',
                6:' -var buff 0.0 -var T 300.0 -var GammaXY 0.2 -var GammaDot 1.0e-04 -var ndump 100 -var ParseData 1 -var DataFile data_init.txt -var DumpFile dumpSheared.xyz',
                4:' -var T 600 -var t_sw 20.0 -var DataFile Equilibrated_600.dat -var nevery 1000 -var ParseData 1 -var WriteData swapped_600.dat', 
                7:' -var buff 0.0 -var T 300 -var nevery 1000 -var ParseData 1 -var DataFile data_minimized.dat -var DumpFile dumpThermalized.xyz -var WriteData equilibrated.dat',
                8:' -var buff 3.0 -var T 0.1 -var sigm 1.5 -var sigmdt 0.01 -var ParseData 1 -var DataFile Equilibrated_300.dat -var DumpFile dumpSheared.xyz',
                9:' -var natoms 1000 -var cutoff 3.52 -var ParseData 1',
                10:' -var T 300.0 -var teq	1.0	-var up -1.0e-03 -var nevery 50 -var ParseData 1 -var DataFile data.0.txt -var DumpFile dumpUp_',
                101:' -var T 300.0 -var teq	1.0	-var up 1.0e-03 -var nevery 50 -var ParseData 1 -var DataFile data.0.txt -var DumpFile dumpDown_',
                11:' -var T 300.0 -var A 0.1 -var Tp 10.0 -var nevery 100 -var DumpFile shearOscillation.xyz -var ParseData 1 -var DataFile data.0.txt', #--- temp(T), amplitude in distance (A), period (Tp)
                5:' -var buff 0.0 -var nevery 1000 -var ParseData 0 -var natoms 1000 -var ntype 2 -var cutoff 3.54  -var DumpFile dumpMin.xyz -var WriteData data_minimized.dat -var seed0 %s -var seed1 %s -var seed2 %s -var seed3 %s'%tuple(np.random.randint(1001,9999,size=4)), 
                13:' -var buff 0.0 -var buffy 0.0 -var T 300 -var swap_every 1 -var swap_atoms 1 -var rn %s -var dump_every 10 -var ParseData 1 -var DataFile equilibrated.dat -var DumpFile traj.dump'%np.random.randint(1001,100000),
                'p0':' swapped_600.dat 10.0 %s'%(os.getcwd()+'/../postprocess'),
                'p1':' swapped_600.dat ElasticConst.txt DumpFileModu.xyz %s'%(os.getcwd()+'/../postprocess'),
                'p2':' %s 3.52 40.0 20.0 40.0 data.txt'%(os.getcwd()+'/../postprocess'),
                } 
    #--- different scripts in a pipeline
    indices = {
                2:[10,101], #--- elastic moduli at finite T
                3:[11], #--- elastic moduli at finite T: laos
                1:[5,7,0,6], #--- minimize, thermalize,melt & quench, shear
                0:[5,7,0,13], #--- minimize, thermalize,melt & quench, anneal
              }[0]
    Pipeline = list(map(lambda x:LmpScript[x],indices))
    Variables = list(map(lambda x:Variable[x], indices))
    EXEC = list(map(lambda x:'lmp' if type(x) == type(0) else 'py', indices))	
    #
    DeleteExistingFolder = True
    #
    EXEC_lmp = ['lmp_mpi','lmp_serial'][0]
    durtn = ['96:59:59','01:59:59'][1]
    mem = '8gb'
    partition = ['gpu-v100','parallel','cpu2019','single'][1]
    #---
    if DeleteExistingFolder:
        os.system( 'rm -rf %s' % jobname ) #--- rm existing
    os.system( 'rm jobID.txt' )
    # --- loop for submitting multiple jobs
    irun = 0
    for counter in runs:
#               cutoff = cutoffs[ irun ]
        print ' i = %s' % counter
        writPath = os.getcwd() + '/%s/Run%s' % ( jobname, counter ) # --- curr. dir
        os.system( 'mkdir -p %s' % ( writPath ) ) # --- create folder
        if irun == 0: #--- cp to directory
            path=os.getcwd() + '/%s' % ( jobname)
            os.system( 'ln -s %s/%s %s' % ( EXEC_DIR, EXEC_lmp, path ) ) # --- create folder & mv oar scrip & cp executable
        #---
        for script,indx in zip(Pipeline,range(100)):
#			os.system( 'cp %s/%s %s/lmpScript%s.txt' %( SCRPT_DIR, script, writPath, indx) ) #--- lammps script: periodic x, pxx, vy, load
            os.system( 'ln -s %s/%s %s' %( SCRPT_DIR, script, writPath) ) #--- lammps script: periodic x, pxx, vy, load
        if sourceFiles: 
            for sf in sourceFiles:
                os.system( 'ln -s %s/Run%s/%s %s' %(sourcePath, counter, sf, writPath) ) #--- lammps script: periodic x, pxx, vy, load
        #---
        makeOAR( path, 1, nThreads, durtn) # --- make oar script
        os.system( 'chmod +x oarScript.sh; mv oarScript.sh %s' % ( writPath) ) # --- create folder & mv oar scrip & cp executable
        os.system( 'sbatch --partition=%s --mem=%s --time=%s --job-name %s.%s --output %s.%s.out --error %s.%s.err \
                            --chdir %s -c %s -n %s %s/oarScript.sh >> jobID.txt'\
                           % ( partition, mem, durtn, jobname[:4], counter, jobname[:4], counter, jobname[:4], counter \
                               , writPath, nThreads, nNode, writPath ) ) # --- runs oarScript.sh! 
        irun += 1

    os.system( 'mv jobID.txt %s' % ( os.getcwd() + '/%s' % ( jobname ) ) )
