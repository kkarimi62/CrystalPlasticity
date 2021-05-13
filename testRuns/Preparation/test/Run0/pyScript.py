import imp
gn=imp.load_source('generate.name','/home/kamran.karimi1/Project/git/CrystalPlasticity/py/generate.py')
gn.Generate( 50688, 5, 0.1,title = 'data.txt',ratio1 = 0.05, ratio2 = 0.26, ratio3 = 0.02, ratio4 = 0.4, ratio5 = 0.27 )
