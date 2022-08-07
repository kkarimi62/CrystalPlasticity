#--- system libraries
import pdb
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
import traceback
import os
import scipy.interpolate as scp_int
import warnings
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import patches
import sys
import sklearn
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_validate
import patsy
import sklearn
from sklearn import linear_model, mixture
import sklearn.mixture as skm
from scipy import optimize
import scipy
import re
from scipy.stats import norm, kurtosis
from pandas.plotting import scatter_matrix
from scipy.ndimage import gaussian_filter
import time
from scipy import ndimage
from scipy.stats import chi2
from scipy.optimize import curve_fit
#--- user-defined
import LammpsPostProcess2nd as lp
import imp
imp.reload(lp)

def GetQuantile(df,q):
    s=df.to_list()
    s.sort()
    n=len(s)
    return s[int(q*n)]

def Wrapper_neighList(lmpData,reference_frames,cutoff):
    '''
    fetch neighbor list from ovito
    '''
    
    fileRef = 'neighList/dump_ref.xyz'
    output  = 'neighList/neighList.xyz'
    #--- rm existig file
    os.system('rm %s'%output)
    #--- split dump file
    for ii0 in reference_frames:
        atom_reference = lp.Atoms(**lmpData.coord_atoms_broken[ii0])
        box0 = lp.Box( BoxBounds = lmpData.BoxBounds[ii0], AddMissing = np.array([0.0,0.0,0.0] ))
        lp.WriteDumpFile(atom_reference, box0).Write(fileRef, itime=ii0,
                 attrs=['id', 'type','x', 'y', 'z'],
                 fmt='%i %i %15.14e %15.14e %15.14e')
        #--- load to ovito
        os.system('ovitos OvitosCna.py %s %s 1 4 %s'%(fileRef,output,cutoff))
        #--- concat
        
        
def WrapperD2min(lmpData,reference_frames,current_frames, dim=3):
    '''
    invoke d2min analysis in ovito
    '''
    #--- split dump file
    for ii0, ii in zip(reference_frames,current_frames):
        atom_current = lp.Atoms(**lmpData.coord_atoms_broken[ii])
        atom_reference = lp.Atoms(**lmpData.coord_atoms_broken[ii0])
        box  = lp.Box( BoxBounds = lmpData.BoxBounds[ii],  AddMissing = np.array([0.0,0.0,0.0] ))
        box0 = lp.Box( BoxBounds = lmpData.BoxBounds[ii0], AddMissing = np.array([0.0,0.0,0.0] ))
        if dim == 3:
        	lp.WriteDumpFile(atom_current, box).Write('D2minAnl/dump_curr.xyz', itime = ii,
                 attrs=['id', 'type','x', 'y', 'z'],
                 fmt='%i %i %15.14e %15.14e %15.14e')
        	lp.WriteDumpFile(atom_reference, box0).Write('D2minAnl/dump_ref.xyz', itime=ii0,
                 attrs=['id', 'type','x', 'y', 'z'],
                 fmt='%i %i %15.14e %15.14e %15.14e')
        else:
        	lp.WriteDumpFile(atom_current, box).Write('D2minAnl/dump_curr.xyz', itime = ii,
                 attrs=['id', 'type','x', 'y'],
                 fmt='%i %i %15.14e %15.14e')
        	lp.WriteDumpFile(atom_reference, box0).Write('D2minAnl/dump_ref.xyz', itime=ii0,
                 attrs=['id', 'type','x', 'y'],
                 fmt='%i %i %15.14e %15.14e')
    #    os.system('tar czf dump.gz dump.xyz')
        fileCurr = 'D2minAnl/dump_curr.xyz'
        fileRef = 'D2minAnl/dump_ref.xyz'
        output = 'D2minAnl/d2min.%s.xyz'%ii
        #--- load to ovito
        os.system('ovitos OvitosCna.py %s %s 2 2 %s'%(fileCurr,output,fileRef))
        #--- concat
        os.system('cat %s >> D2minAnl/d2min.xyz;rm %s'%(output,output))

def WrapperStrain(lmpData,reference_frames,current_frames, dim=3):
    '''
    invoke strain analysis in ovito
    '''
    #--- split dump file
    for ii0, ii in zip(reference_frames,current_frames):
        atom_current = lp.Atoms(**lmpData.coord_atoms_broken[ii])
        atom_reference = lp.Atoms(**lmpData.coord_atoms_broken[ii0])
        box  = lp.Box( BoxBounds = lmpData.BoxBounds[ii],  AddMissing = np.array([0.0,0.0,0.0] ))
        box0 = lp.Box( BoxBounds = lmpData.BoxBounds[ii0], AddMissing = np.array([0.0,0.0,0.0] ))
        if dim == 3:
            lp.WriteDumpFile(atom_current, box).Write('strain/dump_curr.xyz', itime = ii,
                     attrs=['id', 'type','x', 'y', 'z'],
                     fmt='%i %i %15.14e %15.14e %15.14e')

            lp.WriteDumpFile(atom_reference, box0).Write('strain/dump_ref.xyz', itime=ii0,
                 attrs=['id', 'type','x', 'y', 'z'],
                 fmt='%i %i %15.14e %15.14e %15.14e')
        else:
            lp.WriteDumpFile(atom_current, box).Write('strain/dump_curr.xyz', itime = ii,
                     attrs=['id', 'type','x', 'y'],
                     fmt='%i %i %15.14e %15.14e')

            lp.WriteDumpFile(atom_reference, box0).Write('strain/dump_ref.xyz', itime=ii0,
                 attrs=['id', 'type','x', 'y'],
                 fmt='%i %i %15.14e %15.14e')
    #    os.system('tar czf dump.gz dump.xyz')
        fileCurr = 'strain/dump_curr.xyz'
        fileRef = 'strain/dump_ref.xyz'
        output = 'strain/d2min.%s.xyz'%ii
        #--- load to ovito
        os.system('ovitos OvitosCna.py %s %s 2 5 %s'%(fileCurr,output,fileRef))
        #--- concat
        os.system('cat %s >> strain/strain.xyz;rm %s'%(output,output))
        
def WrapperDisp(lmpData,reference_frames,current_frames, dim=3):
    '''
    invoke strain analysis in ovito
    '''
    #--- split dump file
    for ii0, ii in zip(reference_frames,current_frames):
        atom_current = lp.Atoms(**lmpData.coord_atoms_broken[ii])
        atom_reference = lp.Atoms(**lmpData.coord_atoms_broken[ii0])
        box  = lp.Box( BoxBounds = lmpData.BoxBounds[ii],  AddMissing = np.array([0.0,0.0,0.0] ))
        box0 = lp.Box( BoxBounds = lmpData.BoxBounds[ii0], AddMissing = np.array([0.0,0.0,0.0] ))
        lp.WriteDumpFile(atom_current, box).Write('disp/dump_curr.xyz', itime = ii,
                     attrs=['id', 'type','x', 'y', 'z'],
                     fmt='%i %i %15.14e %15.14e %15.14e')

        lp.WriteDumpFile(atom_reference, box0).Write('disp/dump_ref.xyz', itime=ii0,
                 attrs=['id', 'type','x', 'y', 'z'],
                 fmt='%i %i %15.14e %15.14e %15.14e')
        fileCurr = 'disp/dump_curr.xyz'
        fileRef = 'disp/dump_ref.xyz'
        output = 'disp/d2min.%s.xyz'%ii
        #--- load to ovito
        os.system('ovitos OvitosCna.py %s %s 2 6 %s'%(fileCurr,output,fileRef))
        #--- concat
        os.system('cat %s >> disp/disp.xyz;rm %s'%(output,output))
        
def GetAtoms( filee, nevery = 1 ):
    lmpData = lp.ReadDumpFile( filee )
    lmpData.GetCords( ncount = sys.maxsize, 
                     columns = {'c_peratom[1]':'sxx','c_peratom[2]':'syy','c_peratom[3]':'szz',
                                                       'c_peratom[4]':'sxy','c_peratom[5]':'sxz','c_peratom[6]':'syz'}
                    )
    nn=len(lmpData.coord_atoms_broken.keys())
    itimee=list(lmpData.coord_atoms_broken.keys())[0:nn:nevery]

    #--- volume
    rad1=0.0#AtomicRadius[1]
    rad2=0.0#AtomicRadius[2]
    rad3=0.0#AtomicRadius[3]
    neveryy = 10*nn #--- only initial frame for computing volumes
    os.system('ovitos OvitosCna.py %s %s %s %s %s %s %s'%(filee,'Voronoi.xyz',neveryy,3,rad1,rad2,rad3))  
#--- read from d2min.xyz
    ovtData = lp.ReadDumpFile( 'Voronoi.xyz' )
    ovtData.GetCords( ncount = sys.maxsize)
    #--- atom obj
    box0 = map(lambda x:lp.Box( BoxBounds = lmpData.BoxBounds[x], AddMissing = np.array([0.0,0.0,0.0] ) ), itimee ) #--- reference state
    atoms0 = map(lambda x: lp.Atoms( **lmpData.coord_atoms_broken[x].to_dict(orient='series'),
                        AtomicVolume = ovtData.coord_atoms_broken[0]['AtomicVolume'].tolist()), itimee )

    return dict(zip(itimee,list(atoms0))), dict(zip(itimee,list(box0)))

def GetTimeAverageSeries( atoms, col='x' ):
    sarr = np.array(list(map( lambda x: atoms[x][col], atoms.keys() ))).T
    return np.mean(sarr,axis=1)

def GetTimeAverageAtom( atoms0 ):
#    pdb.set_trace()
    itime0 = list(atoms0.keys())[0]
    attrs=atoms0[itime0].__dict__.keys()
    #
    sarr = np.array(list(map(lambda x: GetTimeAverageSeries( atoms0, col=x ), attrs))).T
    df = pd.DataFrame(sarr,columns=attrs)
    #--- pbc: mean of x is not valid! set it to initial values
#    pdb.set_trace()
    df.x = atoms0[itime0]['x'].to_list()
    df.y = atoms0[itime0]['y'].to_list()
    df.z = atoms0[itime0]['z'].to_list()
    #---
    sdict= df.to_dict(orient='series')    
    atomi=lp.Atoms(**sdict)
#    pd.DataFrame(atomi.__dict__)
    return atomi


def GetFrames( lmpData, **kwargs):
    '''
    return frame index
    '''
    if 'nevery' in kwargs and kwargs['nevery'] > 0: #--- return indices nevery times
        return np.arange(0,len(lmpData.coord_atoms_broken.keys()),kwargs['nevery'])
    if 'times' in kwargs: #--- return frame index
        timee = np.array(list(lmpData.coord_atoms_broken.keys()))
        nframe = len(lmpData.coord_atoms_broken.keys())
        return np.array(list(map(lambda x:np.arange(nframe)[timee == x][0],kwargs['times'])))
    
    
def PdfD2min( d2min, 
             #times,     
             Plott = True,
             **kwargs
             ):
    #--- plot
    if Plott:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        #    ax.set_yscale('log')
        #ax.set_xscale('log')
        # ax.set_ylim(1e-5,10)
        # ax.set_xlim(1e-2,1e4)
        ax.set_xlabel(r'log$D^2$min($A^2$)',fontsize=16)
        ax.set_ylabel(r'PDF',fontsize=16)
        ax.tick_params(labelsize=16)
        #ax.set_title(r'itime=%s'%itime)

#     Mean = []
#     Std = []
#     Ebulk = []
#     D2min = {}
    if 1:
#    for itimee in sorted(times): #[0::nn]:

        # strain
#        ebulk = GetStrain(lmpData, [itimee], 0 )[itimee]
#         if ebulk == 0.0:
#             continue
#        Ebulk += [ebulk]    

        #--- d2min
#        d2min = lp.Atoms( **lmpData.coord_atoms_broken[itimee].to_dict(orient='list') )
#        D2min[ itimee ] = np.array(d2min.d2min)
#        Std += [ np.std(d2min.d2min) ]
#        Mean += [ np.mean(d2min.d2min) ]

    #--- size distribution
        if Plott:
            value = np.log10(d2min.d2min)
            #--- filter based on the given limit
#             if 'Limit' in kwargs:
#                 (xlo,xhi)=kwargs['Limit']
#                 filtr = np.all([value>=xlo,value<xhi],axis=0)
#                 value = value[filtr]
#                bins=np.linspace(xlo,xhi,)
            #--- histogram
            hist, edges2, error = GetPDF( value, linscale = True, n_per_decade=32)

        #
            ax.errorbar(edges2,hist,error,fmt='-o',
                        markersize=8,markeredgewidth=0.7,
                            linewidth=.5,
                             barsabove=None,capsize=5,capthick=1,elinewidth=1) #,label='%3.2f'%ebulk)
    #
    if Plott:
#        ax.legend(frameon=False, fontsize=12)
        plt.savefig('pdfD2min.png',dpi=75,bbox_inches='tight')
        plt.show()
        
 #   return Ebulk, Mean, Std, D2min

#def FilterDataFrame(df,column,limits):
#    (xlo,xhi) = limits
#    filtr = np.all([df[column]>=xlo,df[column]<xhi],axis=0)
#    return df[filtr]


def PlotNonLinearDecisionBoundary( ax, X, clf,  ngrid  ):
    #--- draw decison boundary
    xx=np.linspace(np.min(X[:,0]),np.max(X[:,0]),ngrid) #--- grid
    yy=np.linspace(np.min(X[:,1]),np.max(X[:,1]),ngrid)
    xv, yv = np.meshgrid(xx, yy)
    n_grid = xv.shape[0] * xv.shape[1]
#    pdb.set_trace()
    #--- response on grid points 
    XX=np.c_[xv.flatten(),yv.flatten()]
#    hmat = poly.fit_transform(XX)

    z = clf.predict( XX )
    
    #--- reshape z
    z = z.reshape((xv.shape[0],xv.shape[1]))
    
    
    #--- plot

    ax.contour(xv,yv,z,[0.0],colors='red')
    
    plt.show()
    
    
def FetchData(ax, item ):
    xdata = np.array(ax.lines[item].get_xdata())
    ydata = np.array(ax.lines[item].get_ydata())
    yerr = ydata - np.array(ax.lines[item+1].get_ydata())
    #--- sort
    slist = list(zip(xdata,ydata,yerr))
    slist.sort()

    xdata = np.array(slist)[:,0]
    ydata = np.array(slist)[:,1]
    yerr = np.array(slist)[:,2]

    return xdata, ydata, yerr

def Filtr(xdata,xlo,xhi):
    return np.all([xdata>xlo,xdata<xhi],axis=0)


class Stats:
    #--------------------------
    #--- cluster statistics
    #--------------------------
    def __init__( self, mask, xlin, ylin, zlin,
              verbose = False ):
        self.mask = mask
        self.xlin = xlin
        self.ylin = ylin
        self.zlin = zlin
#        if verbose:
#            print('p=%s\npinf=%s\nsmean=%s\nsi_sq=%s'%(p,pinf,smean,crltnl_sq))
        
    def GetProbInf(self):
#        percCount = self.stats[self.stats['percTrue']==True].shape[0]
#        self.pinf0 = 1.0*percCount/self.stats.shape[0]
        
        self.pinf = (self.stats['percTrue'] * self.stats['size']).sum()/self.stats['size'].sum()
#        print(self.pinf0,self.pinf)
    #--- p
    def GetProb(self):
        (ny,nx,nz) = self.mask.shape
        nsize = nx*ny*nz
        self.p = 1.0*self.mask.sum()/nsize #--- occupation prob.

    #--- <s^2>/<s>
    def GetSmean(self):
        self.smean = (self.stats['size']*self.stats['size']).sum()/self.stats['size'].sum()
    #--- correlation length
    def GetCrltnLenSq(self):
        self.si_sq = 2*(self.stats['rg_sq']*self.stats['size'] * self.stats['size']).sum()/\
                  (self.stats['size'] * self.stats['size']).sum()       

    def isPercolating(self,sliceX,sliceY,sliceZ,size):
        (ny,nx,nz)=size
        #
        xlo = sliceX.start
        xhi = sliceX.stop
        assert xhi - xlo <= nx
        #    
        ylo = sliceY.start
        yhi = sliceY.stop
        assert yhi - ylo <= ny
        #    
        zlo = sliceZ.start
        zhi = sliceZ.stop
        assert zhi - zlo <= nz
        #
        return xhi - xlo == nx or yhi - ylo == ny or zhi - zlo == nz

    def GetSize(self):
        #--- clusters
        label_im, nb_labels = ndimage.label(self.mask)
        self.label_im = label_im
        #--- cluster bounds
        sliced=ndimage.find_objects(label_im,max_label=0)
    #    sliceX = sliced[0][1]
    #    sliceY = sliced[0][0]
    #    sliceZ = sliced[0][2]
    #    isPercolating(sliceX,sliceY,sliceZ,mask.shape)
    #
        #--- percolation
        percTrue = list(map(lambda x:self.isPercolating(x[1],x[0],x[2],self.mask.shape),sliced))    
        assert len( percTrue ) == nb_labels

        #--- geometry
        xc = ndimage.measurements.center_of_mass(self.mask, label_im,np.arange(1, nb_labels+1)) #--- (yc,xc,zc)

        (ny,nx,nz) = self.mask.shape
        xv,yv,zv=np.meshgrid(range(nx),range(ny),range(nz))
        xc=ndimage.mean(xv, label_im, np.arange(1, nb_labels+1))
        yc=ndimage.mean(yv, label_im, np.arange(1, nb_labels+1))
        zc=ndimage.mean(zv, label_im, np.arange(1, nb_labels+1))
        varx=ndimage.variance(xv, label_im, np.arange(1, nb_labels+1))
        vary=ndimage.variance(yv, label_im, np.arange(1, nb_labels+1))
        varz=ndimage.variance(zv, label_im, np.arange(1, nb_labels+1))
        #---
        dx = self.xlin[1]-self.xlin[0]
        dy = self.ylin[1]-self.ylin[0]
        dz = self.zlin[1]-self.zlin[0]
        radg_sq = varx * dx * dx + vary * dy * dy + varz * dz * dz
        #
        ones = np.ones(nx*ny*nz).reshape(ny,nx,nz)
        size=ndimage.sum(ones, label_im, np.arange(1, nb_labels+1)) * dx * dy * dz


        #--- postprocess
        df=pd.DataFrame(np.c_[label_im.flatten()],columns=['id'])
        sdict=df.groupby(by='id').groups

    #    pdb.set_trace()
        df_cls = pd.DataFrame(np.c_[list(sdict.keys()),
    #                                list(map(lambda x:len(sdict[x]),sdict.keys()))
                                   ],columns=['cls_id'])

        #--- sort based on id
        df_cls.sort_values('cls_id',ascending=True,inplace=True)
        df_cls = df_cls.iloc[1:] #--- remove label 0

        #--- append percTrue
#        df_cls=pd.DataFrame(np.concatenate((np.c_[df_cls],np.c_[size],np.c_[radg_sq],np.c_[percTrue]),axis=1, dtype=np.object), columns=['cls_id','size','rg_sq','percTrue'])
        df_cls=pd.DataFrame(np.concatenate((np.c_[df_cls],np.c_[size],np.c_[radg_sq],np.c_[percTrue]),axis=1), columns=['cls_id','size','rg_sq','percTrue'])
#        pdb.set_trace()
        #---
        #--- sort based on size
        df_cls.sort_values('size',ascending=False,inplace=True)


        self.stats = df_cls

    def Print(self,xpath,filee,x, header):
            os.system('mkdir -p %s'%xpath)
            np.savetxt('%s/%s'%(xpath,filee),x,header=header)
            

def funcc1(x,a,alpha):
     return a*x**alpha

def funcc0(x,a,beta, alpha,rc):
    return a/((x/rc)**beta+(x/rc)**alpha)

def funcc2(x,a, alpha,rc):
    return a/(1.0+(x/rc)**alpha)
def funcc3(x,a, alpha,rc):
    return a*np.exp(-(x/rc)**alpha)



def nonlinfit(x,y,func_indx,**kwargs):

    #--- initialize 
    if 'yerr' in kwargs: 
        ye=kwargs['yerr']
    #--- sort
    ss=list(zip(x,range(len(x))))
    ss.sort()
    indx=np.array([i[1] for i in ss])
    x=x[indx]
    y=y[indx]
    if 'yerr' in kwargs: 
        ye = ye[indx]
    #--- filter:
    y=y[x>0]
    if 'yerr' in kwargs: 
        ye=ye[x>0]
    x=x[x>0]
    #
    #--- fit
    if func_indx == 0:
        func = funcc0
        p0=(1.0e4,1.0,1.12,0.1) #a,beta, alpha,rc)
        bounds=((0.0,0.0,0.0,0.0), (np.infty,np.infty,np.infty,0.5))   
    elif func_indx == 1:
        func = funcc1
        p0=(1.0e4,1.12) #(a,alpha)
        bounds=((0.0,0.0), (np.infty,np.infty))
    elif func_indx == 2:
        func = funcc2
        p0=(1.0e4,1.12,0.1) #(a,alpha,rc)
        bounds=((0.0,0.0,0.0), (np.infty,np.infty,np.infty))
        if 'rc' in kwargs:
            rc_range=kwargs['rc']
            p0=(1.0e4,1.12,np.mean(rc_range)) #(a,alpha,rc)
            bounds=((0.0,0.0,rc_range[0]), (np.infty,np.infty,rc_range[1]))
    elif func_indx == 3:
        func = funcc3
        p0=(1.0e4,1.12,0.1) #(a,alpha,rc)
        bounds=((0.0,0.0,0.0), (np.infty,np.infty,np.infty))
        print('hello!')

    if 'bounds' in kwargs:
        bounds = kwargs['bounds']
    if 'p0' in kwargs:
        p0 = kwargs['p0']
        
    if func_indx == 1 and 'fit' in kwargs and kwargs['fit'] == 'linear' :
#        assert not np.any(np.log10(y)==-np.infty), 'y<=0!'
        #--- filter
        ye=ye[y>0]
        x=x[y>0]
        y=y[y>0]

        popt, pcov = np.polyfit(np.log10(x),np.log10(y),deg=1, w=1.0/ye,cov=True)
        alpha = popt[0]
        a = 10**popt[1]
        popt=[a,alpha]
        
    else:
        popt, pcov = curve_fit(func, x, y, 
                           p0=p0,
                           bounds=bounds,    
                           sigma=ye)

    # #--- goodness of fit
    # n=len(x)
    # m=4 #--- number of coefficients
    # xsq = np.sum(((y-func(x,*popt0))/ye)**2)
    # pval=1-chi2.cdf(xsq,n-m)

    # if pval > pval_max:
    #     pval_max = pval
    
    return popt, pcov


def exponents(xdata,ydata, verbose=False, Plot = False, **kwargs):
    if len(ydata[~np.isnan(ydata)]) == 0:
        return np.nan, np.nan, np.nan

    if Plot:
        #--- plot scatter
        ax = PltErr(xdata,ydata, 
           Plot = False,
            attrs={'fmt':'.', 'markersize':'12'},
           **kwargs
              )

    #--- compute rc: the fit is a/((x/rc)**beta+(x/rc)**alpha)
    # xdata = np.array(ax_s.lines[0].get_xdata())
    # ydata = np.array(ax_s.lines[0].get_ydata())
    # popt0=\
    # nonlinfit( xdata,
    #            ydata,
    #           0,
    # #           yerr = None
    #          )
    # (a,beta,alpha,rc) = popt0

    #--- compute rc, fit: a/(1.0+(x/rc)**alpha)
    if not 'rc' in kwargs:
        popt2,pcov=\
        nonlinfit( xdata,
                   ydata,
                  2, #--- fit index
                  **kwargs,
                 )
        (a,alpha,rc) = popt2
        
        if Plot:


            #--- plot fit
            ax=PltErr(xdata,funcc2(xdata,*popt2), 
               Plot = False,
                attrs={'fmt':'-.',
                       'color':'red',
                        'label':r'$x^{-%3.2f},r_c=%2.1e$'%(alpha,rc),

                      },
                ax=ax,
                    legend=True,
                  )
        
        
    else:
        rc = kwargs['rc']
    
    if Plot:
         #--- plot rc
        ax=PltErr([rc,rc],ax.axis()[2:], #[ax.axis()[2],funcc2(rc,*popt2)], 
           Plot = False,
            attrs={'fmt':'-.',
                   'color':'red',

                  },
            ax=ax,
              )
    #--- max curvature
#     xxdata=xdata[xdata>0.0]
#     Xdata=np.logspace(np.log10(xxdata.min()),np.log10(xxdata.max()),128)
#     Ydata = funcc2(Xdata,*popt2)
#     curv = np.gradient(np.gradient(np.log10(Ydata),np.log10(Xdata)),np.log10(Xdata))
#     slist = [[i,j] for i, j in zip(curv,Xdata)]
#     slist.sort(reverse=False)
#     print(curv)
#     rc = slist[0][1]



#    if verbose: 
#        print('rc=',rc,'alpha=',alpha)
        
        
        


    #--- fit a/x**alpha within the range x>rc
    xdataf = xdata[xdata>rc]
    ydataf = ydata[xdata>rc]
    yerr = kwargs['yerr'][xdata>rc]
    
    if verbose:
        print('xdataf=',xdataf)

    
    popt1,pcov=\
    nonlinfit( xdataf,
               ydataf,
              1,
              fit = 'linear',
              yerr=yerr,
             )
    (a,alpha) = popt1
    if verbose: 
        print('alpha=',alpha)
        print('err=',pcov[0,0]**0.5)
        
    if Plot and 'rc' in kwargs:


        #--- plot fit
        PltErr(xdata,funcc1(xdata,*popt1), 
           Plot = False,
            attrs={'fmt':'-.',
                   'color':'red',
                    'label':r'$x^{%3.2f},r_c=%2.1e$'%(alpha,rc),

                  },
            ax=ax,
                legend=True,
              )        

    return a, rc, alpha, pcov[0,0]**0.5


def exponents2nd(xdata,ydata, verbose=False, Plot = False, **kwargs):
    if len(ydata[~np.isnan(ydata)]) == 0:
        return np.nan, np.nan, np.nan
    #--- compute rc: the fit is a/((x/rc)**beta+(x/rc)**alpha)
    # xdata = np.array(ax_s.lines[0].get_xdata())
    # ydata = np.array(ax_s.lines[0].get_ydata())
    # popt0=\
    # nonlinfit( xdata,
    #            ydata,
    #           0,
    # #           yerr = None
    #          )
    # (a,beta,alpha,rc) = popt0

    #--- compute rc, fit: a/(1.0+(x/rc)**alpha)
        
#     popt2=\
#     nonlinfit( xdata,
#                ydata,
#               2,
#     #           yerr = None
#               **kwargs,
#              )
#     (a,alpha,rc) = popt2
#     if verbose: 
#         print('rc=',rc,'alpha=',alpha)
        
        
        

    rc = kwargs['rc']
    #--- fit a/x**alpha within the range x>rc
    #rc = 0.1 #!!!!1
    xdataf = xdata[xdata>rc]
    ydataf = ydata[xdata>rc]
    popt1=\
    nonlinfit( xdataf,
               ydataf,
              1,
              fit = 'linear',
    #           yerr = None
             )
    (a,alpha) = popt1
    if verbose: 
        print('alpha=',alpha)
        
    if Plot:
        #--- plot scatter
        ax = PltErr(xdata,ydata, 
           yerr=None, 
           Plot = False,
            attrs={'fmt':'.', 'markersize':'12'},
           **kwargs
              )

        #--- plot fit
        ax=PltErr(xdataf,funcc1(xdataf,*popt1), 
           yerr=None, 
           Plot = False,
            attrs={'fmt':'-.',
                   'color':'red',
                    'label':r'$x^{-%3.2f},r_c=%2.1e$'%(alpha,rc),

                  },
            ax=ax,
                legend=True,
           **kwargs
              )
         #--- plot rc
#         ax=PltErr([rc,rc],[ax.axis()[2],funcc2(rc,*popt2)], 
#            yerr=None, 
#            Plot = False,
#             attrs={'fmt':'-.',
#                    'color':'red',

#                   },
#             ax=ax,
#            **kwargs
#               )

        
    return a, rc, alpha


def Pltt(ax_s,ax_m,**kwargs):
    #--- fetch scatter data
    X = np.c_[ax_s.lines[0].get_xdata(),ax_s.lines[0].get_ydata()]
    if 'indices' in kwargs: #--- concat
        xdata = np.array(list(map(lambda x:ax_s.lines[x].get_xdata(),kwargs['indices']))).flatten()
        ydata = np.array(list(map(lambda x:ax_s.lines[x].get_ydata(),kwargs['indices']))).flatten()
        X = np.c_[xdata,ydata]
#        print(X)

    #--- plot
    attrs={ 'color':'C0',
            'markerfacecolor':'C0',
            'markeredgecolor':'white',
            'markersize':10,
            'marker':'o',
           'markeredgewidth':1.75,
            'linewidth':1, 
           'barsabove':None,
           'capsize':5,
           'capthick':1,
           'elinewidth':1,
           'fmt':'.',
#           'markevery':2,
          }       
    ax = None  
    ax = PltErr(X[:,0],X[:,1], 
       yerr=None, 
       Plot = False,
        attrs = attrs,
#        zorder=0,
      )




    #--- plot average
    Xm = np.c_[ax_m.lines[0].get_xdata(),
               ax_m.lines[1].get_ydata(),
               np.array(ax_m.lines[2].get_ydata())-np.array(ax_m.lines[0].get_ydata())]


    attrs={ 'color':'red',
            'markersize':10,
            'marker':'s',
            'markerfacecolor':'red',
            'markeredgecolor':'white',
           'markeredgewidth':1.75,
            'linewidth':1, 
           'barsabove':None,
           'capsize':5,
           'capthick':1,
           'elinewidth':1,
          }  
    ax=PltErr(Xm[:,0],Xm[:,1], 
            yerr=Xm[:,2], 
            Plot = False,
#            DrawFrame=DFset, #(0.17,0.17,0.15,0.06,0.1),
            attrs = attrs,
            ax=ax,
#            zorder=5,
            **kwargs,
#           yticks=(['0','10','20'],['0','10000','20000'])
      )


    
    return ax






def GetMismatchIco(atoms,box,
                   df_comp,
                   xlin,ylin,zlin,
                  ):
    nx,ny,nz = len(xlin), len(ylin), len(zlin)
    #--- partition box & assign index to each atom
    wrap = lp.Wrap(atoms,box)
    wrap.WrapCoord() #--- wrap inside
    wrap.Set(atoms)
    assert np.sum(wrap.isInside()) == len(atoms.x)
    wrap.GetDimensionlessCords()
    AtomCellId = (wrap.beta * np.array([nx,ny,nz])).astype(int)
    #--- store in a df
#    pdb.set_trace()
    tmp=df_comp.set_index(['indxx','indxy','indxz'])
    return tmp.loc[list(AtomCellId)]['deltaa'].to_list()

def GetComp( atoms, atomf ):
    #--- different types
    types = list(set(atomf.type))
    types.sort()
    ntype=len(types)
    c=np.zeros(ntype) #{}
    n = len(atoms.x)
    for typei,indxx in zip(types,range(ntype)):
#        c[typei] = 1.0*np.sum(atoms.type == typei)/n
        c[indxx] = 1.0*np.sum(atoms.type == typei)/n

    return c

def GetDelta(df, AtomicRadius):
    size = list(map(AtomicRadius.get,list(map(int,df['type'].tolist()))))
    assert len(size) == len(df)
    df = pd.DataFrame(np.c_[df,size],columns=list(df.keys())+['size'])
    assert df['size'].mean() > 0.0
    assert df['size'].std() >= 0.0, 'increase size!'
    return df['size'].std()/df['size'].mean()

def GetPressComp( atoms,box, dmean, AtomicRadius, **kwargs ):
#    pdb.set_trace()
    #--- grid: tiling mapped box with original size
    (xlin, ylin, zlin), (xv, yv, zv) = lp.GetCubicGrid( box.CellOrigin, 
                                                     box.CellVector, 
                                                     dmean,
                                                     margin = 0.0 * dmean, odd = False )
    xi = np.array(list(zip(xv.flatten(), yv.flatten(), zv.flatten())))
    dvol = (xlin[1]-xlin[0])*(ylin[1]-ylin[0])*(zlin[1]-zlin[0])
    (ny,nx,nz) = xv.shape
#     nx -= 1
#     ny -= 1
#     nz -= 1
    assert nx*ny*nz >= 8, 'decrease division length!'
    print(dmean,nx*ny*nz)
    #--- indices
    (xvi, yvi, zvi) = np.meshgrid(np.arange(0,nx),np.arange(0,ny),np.arange(0,nz))
    indices = np.array(list(zip(xvi.flatten(), yvi.flatten(), zvi.flatten()))) #--- shape:(ncel,3)
    indices = list(map(lambda x: tuple(x),indices))

    #--- partition box & assign index to each atom
    wrap = lp.Wrap(atoms,box)
    wrap.WrapCoord() #--- wrap inside
    wrap.Set(atoms)
    assert np.sum(wrap.isInside()) == len(atoms.x)
    wrap.GetDimensionlessCords()
    AtomCellId = (wrap.beta * np.array([nx,ny,nz])).astype(int)
    #--- store in a df
    df = pd.DataFrame(np.c_[pd.DataFrame(atoms.__dict__),AtomCellId],
                         columns=list(pd.DataFrame(atoms.__dict__).keys())+['ix','iy','iz'])
    df['ix']=df['ix'].astype(int)
    df['iy']=df['iy'].astype(int)
    df['iz']=df['iz'].astype(int)
#    display(df.head())

    #--- group & compute p and c
    d = df.groupby(by=['ix','iy','iz']).groups
    assert len(d) == nx*ny*nz, 'empty boxes!'
#     #--- lambda function: compute p 
#     f = lambda x: np.sum(np.sum(np.array(x.sxx)+np.array(x.syy)+np.array(x.szz)))*(-1.0e-4/3.0/dvol)
#     vol=np.linalg.det(box.CellVector)
#    pdb.set_trace()
    #
    
    keys = d.keys()
    plist = list(map(lambda x:GetDelta(df.iloc[d[x]],AtomicRadius),keys)) #--- len(plist) = ncell
    clist = list(map(lambda x:GetComp(df.iloc[d[x]],atoms),keys)) #--- clist[icell]={1:c1,2:c2, ...}

    #--- make a data frame
    types = list(set(atoms.type))
    types.sort()
    cols = list(map(lambda x:'Type%s'%x,types))
#    pdb.set_trace()
    df_comp = pd.DataFrame(np.concatenate((np.c_[indices],np.c_[xi],np.c_[clist],np.c_[plist]),
                                          axis=1,dtype=np.object),
                           columns=['indxx','indxy','indxz','x','y','z']+cols+['deltaa'])
    
#    pdb.set_trace()
    #---
    
    
    #--- read file: elastic constants
#     if 'MODU' in kwargs and kwargs['MODU']:
#         fileName = kwargs['PATH']
#         modu = pd.read_csv(fileName, sep=' ',header=0)
# #        display(modu.head())
#         if 'PLOT' in kwargs and kwargs['PLOT']:
#             plott(modu['C66'],nx,ny,nz,box,zlin, 'muCG.png')

    #--- plot
    #--- reshape value


        
#        display(modu.head())
#     if 'MODU' in kwargs and kwargs['MODU']:
#         mlist = modu['C66'].to_list()
#         return clist, plist, mlist
#     else:
    return (xlin,ylin,zlin), df_comp

def PltBinary(xlin,ylin,zlin, 
              val,
              box0,
              thresh = 0.0,
              **kwargs,
             ):
    #--- reshape value
    (nx,ny,nz) = len(xlin), len(ylin),len(zlin) 
    value = np.c_[val].reshape(((ny,nx,nz)))

    #--- mean & variance

    #--- xy plane
    #--- 2d slice
    nzl=[0]
    value2d = Get2dSlice( value, zlin, 
                        zlin[-1], nzll=nzl  )

    #--- square bitmap
    lx=np.min([xlin[-1]-xlin[0],ylin[-1]-ylin[0]])
    xc = 0.5*(xlin[-1]+xlin[0])
    yc = 0.5*(ylin[-1]+ylin[0])                  
    
    #--- plot
    CellVectorOrtho, VectorNorm = lp.GetOrthogonalBasis( box0.CellVector ) #--- box length
    PltBitmap(value2d<thresh, 
#              xlim=VectorNorm[0]*np.array([-0.5,0.5]),ylim=VectorNorm[1]*np.array([-0.5,0.5]),
              xlim=np.array([xc-0.5*lx,xc+0.5*lx]),ylim=np.array([yc-0.5*lx,yc+0.5*lx]),
              frac = 1.0, #--- plot a patch
              **kwargs
            )
    return value

def GetVolume(lmpData,itimee):
    box = lp.Box( BoxBounds = lmpData.BoxBounds[itimee], AddMissing = np.array([0.0,0.0,0.0] ))
    CellVectorOrtho, VectorNorm = lp.GetOrthogonalBasis( box.CellVector )
    return np.linalg.det( CellVectorOrtho )



def row_histogram(a):
    ca = np.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = np.unique(ca, return_index=True, return_inverse=True)
    counts = np.bincount(inverse)
    sort_indices = np.argsort(counts)[::-1]
    return (a[indices[sort_indices]], counts[sort_indices])


        
def GetMotifs( lmpData, times, Plot = True ):
    Unique_indices={}
    Counts = {}
    for itime in times:
        vor = lp.Atoms( **lmpData.coord_atoms_broken[itime].to_dict(orient='list') )
        #
        voro_indices = np.c_[vor.VoronoiIndex3,vor.VoronoiIndex4,vor.VoronoiIndex5,vor.VoronoiIndex6].astype(int)
        # Compute frequency histogram.

        unique_indices, counts = row_histogram(voro_indices)
        Unique_indices[itime]=unique_indices[:]
        Counts[itime]=100.0*counts/len(voro_indices)
        #
        # Print the ten most frequent histogram entries.
        if Plot:
            for i in range(10):
                print("%s\t%i\t(%.1f %%)" % (tuple(unique_indices[i]), 
                                         counts[i], 
                                         100.0*float(counts[i])/len(voro_indices)))
                plt.bar(str(tuple(unique_indices[i])), 100.0*float(counts[i])/len(voro_indices),color='C0')
            #
#             plt.yscale('log')
#             plt.ylim(.1,100)
            plt.xlabel('Voronoi Index(n3,n4,n5,n6)',fontsize=16)
            plt.ylabel('Percentage',fontsize=16)
            plt.tick_params(labelrotation=90,labelsize=16)
            #
            plt.savefig('motifs.png',dpi=75,bbox_inches='tight')
            plt.show()
            
    return Unique_indices, Counts


def PdfCondD2min( lmpData, lmpDmin, times,     
                  Plot = True, title='pdfCondD2min.png',
                 axisLabels = True,
                ):
    #--- plot
    if Plot:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        #
        ax.set_yscale('log')
        #ax.set_xscale('log')
        #
        ax.set_ylim(1,1e4)
        ax.set_xlim(-1,3)
        #
        if axisLabels:
            ax.set_xlabel(r'log$D^2$min($A^2$)',fontsize=16)
            ax.set_ylabel(r'PDF',fontsize=16)
        #
        #ax.set_title(r'itime=%s'%itime)
        #
        PutMinorTicks( ax, LOGY = True)
        #
        ax.tick_params(labelsize=20,which='both',axis='both', top=True, right=True)

    valueFiltrd = {}
    valueTot = {}
    for itimee in sorted(times): 

        # strain
        ebulk = GetStrain(lmpData, itimee, 0 )
#         if ebulk == 0.0:
#             continue

        #--- d2min and vor
        d2min = lp.Atoms( **lmpDmin.coord_atoms_broken[itimee].to_dict(orient='list') )
        vor = lp.Atoms( **lmpData.coord_atoms_broken[itimee].to_dict(orient='list') )

        #--- filter
        indices = GetFullIcosahedra( vor )

        #--- histogram: total
        value = valueTot[itimee]=np.array(d2min.d2min)
        hist, edges, error = GetPDF( np.log10(value), 
                                    linscale = True, 
                                    n_per_decade=32,
                                    density = False,
                                   )
        #--- filtered
        value = valueFiltrd[ itimee ] = np.array(d2min.d2min)[indices]
        hist_filtrd, edges_filtrd, error_filtrd = GetPDF( np.log10(value), 
                                                         linscale = True, 
                                                         n_per_decade=32,
                                                         density = False,
                                                        )
        #--- invert
        valinv = np.array(d2min.d2min)[~indices]
        hist_inv, edges_inv, error_inv = GetPDF( np.log10(valinv), 
                                                linscale = True, 
                                                n_per_decade=32,
                                                density = False,
                                               )


        if Plot:
            attr = { 'markersize':10,'markeredgewidth':0.7,
                     'linewidth':.5,'barsabove':None,'capsize':5,
                     'capthick':1,'elinewidth':1}
            #
    #             ax.errorbar(edges,hist,error,fmt='-o',
    #                         markersize=8,markeredgewidth=0.7,
    #                         linewidth=.5,
    #                          barsabove=None,capsize=5,capthick=1,elinewidth=1,label='Total')
                #
            ax.errorbar(edges_filtrd,hist_filtrd,error_filtrd,
                        **attr,
                        fmt='-o', color='black', 
                        markerfacecolor='white', markeredgecolor=None,
                        label='Icosahedra', markevery = int(len(edges_filtrd)/10),
                        errorevery = int(len(edges_filtrd)/10),
                       )
            #
            ax.errorbar(edges_inv,hist_inv,error_inv,
                        **attr,
                        fmt='-s', color='red',
                        markerfacecolor=None, markeredgecolor='black',
                        label='Non Icosahedra',markevery = int(len(edges_inv)/10),
                        errorevery = int(len(edges_inv)/10),
                       )
    #
    if Plot:
        ax.legend(frameon=False, fontsize=12)
        #
        DrawFrame(ax, 0.2,0.09,0.15,0.06,0.04,LOG_Y=True)
        #
        plt.savefig(title,dpi=2*75,bbox_inches='tight',pad_inches=0.0)
        plt.show()

    return valueFiltrd,valueTot

def Rescale(file0_indx,xdata, ydata, erry):
    if file0_indx in [18]:
        ydata *= xdata ** 1.8
        erry *= xdata ** 1.8

def ReadDataa(file0,file0_indx,verbose=False):
    #--- load data
    xdata = np.array([np.nan])
    ydata = np.array([np.nan])
    erry = np.array([np.nan])
    errx = np.array([np.nan])
    
    if os.path.isfile( file0 ):
        sarr = np.loadtxt(file0)
    else:
        if verbose:
                print('no such a file: %s'%file0)
        pass
    try:
        #--- assign x, y
        if file0_indx in [5]:
            xdata = np.array([sarr[0]])
            ydata = np.array([-sarr[1]])
        elif file0_indx in [ 25, 26, 27,29,30 ]:
            xdata = np.array([sarr[1]])
            ydata = np.array([sarr[0]])
            erry = np.array([sarr[2]]) #--- error
            errx = np.array([sarr[3]]) #--- error
            
        elif file0_indx in [6, 8, 10, 13 ]:
            xdata = np.array([sarr[0]])
            ydata = np.array([sarr[1]])
            
        elif file0_indx in [11, 14, 19, 22]:
            xdata = np.array([sarr[0]])
            ydata = np.array([sarr[2]])
            
        elif file0_indx in [15, 17, 20, 21, 23, 24 ]:
            xdata = np.array([sarr[1]])
            ydata = np.array([sarr[2]])
            
        else:
            xdata = sarr[:,0] #--- x
            ydata = sarr[:,1] #--- y
            
        if file0_indx in [9, 12, 16, 18]:
            erry = sarr[:,2] #--- error

        #--- only for shear band width
        if file0_indx in [3]:
            ydata2 = sarr[:,2] #--- y1
            #--- remove nan's
            indices = np.any(~np.isnan(np.c_[ydata,ydata2]),axis=1) #--- true:either x or y is not nan
            ydata = np.array(list(map(lambda x:np.mean(x[~np.isnan(x)]),np.c_[ydata,ydata2][indices]))) #--- take mean                    
            xdata = xdata[indices]
            assert xdata.shape == ydata.shape
            
        if not len(erry) == len(ydata):
            erry = np.ones(len(ydata))*np.nan
        if not len(errx) == len(xdata):
            errx = np.ones(len(xdata))*np.nan
            
    except:
        if verbose: print('error while reading %s'%file0)
        pass
        #return

    return xdata, ydata, erry, errx
            
            
def PlotPaperVersion(pathh_indx,
                     file0_indx,
                     runs = [0],
                     times = range(0,200+1,2),
                     **kwargs):
    
    verbose = True if 'verbose' in kwargs and kwargs['verbose'] == True else False
    timeseries = True if not 'timeseries' in kwargs else kwargs['timeseries']
    drawstyle = 'default' if not 'drawstyle' in kwargs else kwargs['drawstyle']
    xstr = '' if not 'xlabel' in kwargs else kwargs['xlabel']
    ystr = '' if not 'ylabel' in kwargs else kwargs['ylabel']
    nevery = kwargs['nevery'] if 'nevery' in kwargs else 1
    #--- setup symbols
    colors = ['C0','red','green','blue','cyan','brown','grey','magenta','orange','yellow']
    fillstyles=['white',None,'white',None,'white',None,'white',None,'white',None,'white',None,'white',None,'white',None]
    markers=['o','s','D','^','<','>','v']
    markersizes=[10,10,10,12,12,12,10]
    #--- plot
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    #

    
    for mg, marker, color, fillstyle, markersize in list(zip( [ 
                                         'FeNi',
                                           'CoNiFe',
                                           'CoNiCrFe',
                                            'CoCrFeMn',
                                            'CoNiCrFeMn',
                                            'Co5Cr2Fe40Mn27Ni26',
#                                            'Co5Cr2Fe40Mn27Ni262nd',
#                                             'CuZr3'
                                        ],markers, colors, fillstyles, markersizes )):
        if 'glass' in kwargs and kwargs['glass'] != mg:
            continue
        print(mg)
        Xdata = []
        Ydata = []
        #--- loop over realizations 
        for irun in runs:
            xdata_timeseries = []
            ydata_timeseries = []
            erry_timeseries = []
            errx_timeseries = []
            
            #--- loop over time
            for itimee, count in\
            list(zip(times, range(len(times)))): 

                pathh = { 
                          0:'%s/PairCrltnT300/%s/Run%s'%(os.getcwd(),mg,irun),
                          1:'%s/VorAnlT300/%s/Run%s'%(os.getcwd(),mg,irun),
                          2:'%s/D2minAnalysisT300/%s/Run%s'%(os.getcwd(),mg,irun),
                          3:'%s/ElasticityT300/%s/eps2/itime%s/Run%s'%(os.getcwd(),mg,itimee,irun),
                          4:'%s/ElasticityT300/%s/eps2/itime%s/Run%s/ModuAnl'%(os.getcwd(),mg,itimee,irun),
                          5:'%s/Exponents/%s'%(os.getcwd(),mg),
                        }[pathh_indx]
                file0 = {
                          0:'%s/gr.txt'%pathh,
                          1:'%s/icoPercentageWithStrain.txt'%pathh,
                          2:'%s/d2min_gamma.txt'%pathh,
                          3:'%s/rc_d2min.txt'%pathh,
                          4:'%s/crsD2minRhoSro.txt'%pathh,
                          5:'%s/NegativeModulus.txt'%pathh, 
                          6:'%s/YieldDelta.txt'%pathh, 
                          7:'%s/gr_peak_gamma.txt'%pathh, 
                          8:'%s/muClustersize_gamma.txt'%pathh, 
                          9:'%s/pdfMu.txt'%pathh, 
                          10:'%s/mu_mean_std.txt'%pathh, 
                          11:'%s/mu_mean_std.txt'%pathh, 
                          12:'%s/pdfClsSize.txt'%pathh,
                          13:'%s/muClustersize_gamma.txt'%pathh, 
                          14:'%s/muClustersize_gamma.txt'%pathh, 
                          15:'%s/muClustersize_gamma.txt'%pathh, 
                          16:'%s/ps.txt'%pathh, 
                          17:'%s/muClustersize_gamma.txt'%pathh, 
                          18:'%s/ps.txt'%pathh, 
                          19:'%s/pinf.txt'%pathh, 
                          20:'%s/pinf.txt'%pathh, 
                          21:'%s/pinf.txt'%pathh, 
                          22:'%s/crltnl_gamma.txt'%pathh, 
                          23:'%s/crltnl_gamma.txt'%pathh, 
                          24:'%s/crltnl_gamma.txt'%pathh, 
                          25:'%s/hmin_gamma_exp.txt'%pathh, 
                          26:'%s/hmin_nu_exp.txt'%pathh, 
                          27:'%s/hmin_beta_exp.txt'%pathh, 
                          28:'%s/s_rg.txt'%pathh, 
                          29:'%s/hmin_tau_exp.txt'%pathh, 
                          30:'%s/hmin_df_exp.txt'%pathh, 
                        }[file0_indx]
                
                #--- read data
                xdata, ydata, erry, errx = ReadDataa(file0,file0_indx)
                Rescale(file0_indx,xdata, ydata, erry)

                    
                if not timeseries:   #--- plot each time  
                    if verbose:
                        print('xdata=',xdata,'\nydata=',ydata,'\nerry=',erry,'\nerrx=',errx)                        
                    attrs={ 'color':colors[count],
                            'markersize':markersizes[count],
                            'marker':markers[count],
                            'markerfacecolor':colors[count],
                            'markeredgecolor':'white', #'black' if not fillstyles[count] else None,
                            'label':'%s/irun%s/itime%s'%(mg,irun,itimee),
                           'markevery':nevery,
                           'errorevery':nevery,
                           'markeredgewidth':kwargs['markeredgewidth'] if 'markeredgewidth' in kwargs else 1.75,
                            'linewidth':1, 
                           'barsabove':None,
                           'capsize':5,
                           'capthick':1,
                           'elinewidth':1,
                           'fmt':kwargs['fmt'] if 'fmt' in kwargs else '.',
                           'drawstyle':drawstyle,
                          }
                    #--- rescale
                    if 'scaley' in kwargs:
                        ydata /= kwargs['scaley'][count]
                        erry /= kwargs['scaley'][count]
                    
                    if verbose:
                        print('plot itime=%s,irun=%s,len(ax.lines)=%s'%(itimee,irun,len(ax.lines)))
#                     if erry == [np.nan]:
#                         erry = None
                    PltErr(xdata,ydata, 
                       yerr=erry,      
                       xerr=errx, #None if not 'xerr' in kwargs or not kwargs['xerr'] else  errx,      
                       ax = ax,
                       xstr = xstr,
                       ystr = ystr,
                       attrs = attrs,
                       Plot = False,
                       **kwargs,
                      )

                else:
                    assert timeseries #--- concat different times
                    xdata_timeseries += list(xdata)
                    ydata_timeseries += list(ydata)
                    erry_timeseries  += list(erry)
                    errx_timeseries  += list(errx)

                    
            if timeseries: #--- plot timeseries
                xdata_timeseries = np.array(xdata_timeseries)
                ydata_timeseries = np.array(ydata_timeseries)
                erry_timeseries = np.array(erry_timeseries)
                errx_timeseries = np.array(errx_timeseries)
            #--------------------
            #--- rescale data
            #--------------------
                if file0_indx in [ 15, 24 ]:
                    if np.any(~np.isnan(ydata_timeseries)):
                        pmax = np.max(xdata_timeseries)
                        #
                        tmp = np.array(xdata_timeseries)
                        tmp.sort()
                        q=0.95
                        ns=len(tmp)
                        pmax = tmp[int(q*ns)]
                        #
                        if 'pmax' in kwargs:
                            pmax = kwargs['pmax']

                        if verbose:
                            print('pmax=',pmax)
                        xdata_timeseries = 1.0 - np.array( xdata_timeseries ) / pmax 
                elif file0_indx in [21]:
                    if np.any(~np.isnan(ydata_timeseries)):
#                        pc = xdata_timeseries[ydata_timeseries>0.0][0]
                        ind = np.arange(len(ydata_timeseries))[ydata_timeseries>0.0][0]
                        pc = xdata_timeseries[ind-1]
                        if 'pc' in kwargs:
                            pc = kwargs['pc']
                        if verbose:
                            print('pc=',pc)
                        xdata_timeseries = np.array( xdata_timeseries ) / pc - 1.0

        
                #--- concat different realizations
                try:
                    Xdata = np.concatenate((Xdata,xdata_timeseries),axis=0) #--- concat different realizations
                    Ydata = np.concatenate((Ydata,ydata_timeseries),axis=0) 
                except:
                    traceback.print_exc()
                    Xdata = xdata_timeseries.copy()
                    Ydata = ydata_timeseries.copy()
                

                #--- plot
                if 'PlotEvery' in kwargs and kwargs['PlotEvery']:
                #--- graph-related attribute
                    attrs={ 'color':color,
                            'markersize':markersize,
                            'marker':marker,
#                            'markerfacecolor':fillstyle,
                            'markeredgecolor':'white', #'black' if not fillstyle else None,
                            'label':'%s/irun%s'%(mg,irun),
                           'markevery':nevery,
                           'errorevery':nevery,
                           'markeredgewidth':1.75,
                            'linewidth':1, 
                           'barsabove':None,
                           'capsize':5,
                           'capthick':1,
                           'elinewidth':1,
                           'fmt':kwargs['fmt'] if 'fmt' in kwargs else '.',
                           'zorder':1,
                          }                    
                    #--- plot ydata
                    if verbose:
                        print('ydata=',ydata_timeseries)
                        print('plot irun=%s,len(ax.lines)=%s'%(irun,len(ax.lines)))
                    PltErr(xdata_timeseries,ydata_timeseries, 
                           yerr=erry_timeseries,      
                           xerr=errx_timeseries,      
                           ax = ax,
                           xstr = xstr,
                           ystr = ystr,
                           attrs = attrs,
                           Plot = False,
                           **kwargs,
    #                           xerr=yerr0,
                          )

        #
        if 'PlotMean' in kwargs and kwargs['PlotMean']:
            Xdata = Xdata[~np.isnan(Ydata)]
            Ydata = Ydata[~np.isnan(Ydata)]
            #---
#             if file0_indx in [21]: #--- include positive vals
#                 Xdata = Xdata[Ydata>0]
#                 Ydata = Ydata[Ydata>0]
            if verbose:
                print('ydata=',Ydata)
                print('plot len(ax.lines)=%s'%(len(ax.lines)))
            try:
                #--- take average
                nbins=1024 if not 'nbins' in kwargs else kwargs['nbins']
                if file0_indx == 8:
                    nbins=1
                Xbin, Ybin, Yerr = BinData(Xdata,Ydata,nbins=nbins,scale='linear' if not 'scale' in kwargs else kwargs['scale'])
                attrs={ 'color':color,
                        'markersize':markersize,
                        'marker':marker,
                        'markerfacecolor':fillstyle,
                        'markeredgecolor':'black' if not fillstyle else None,
                        'label':'%s'%(mg),
                       'markevery':nevery,
                       'errorevery':nevery,
                       'markeredgewidth':0.7,
                        'linewidth':1, 
                       'barsabove':None,
                       'capsize':5,
                       'capthick':1,
                       'elinewidth':1,
                       'zorder':2,
                      }

                PltErr(Xbin,Ybin, 
                       yerr=Yerr,      
                       ax = ax,
                       attrs = attrs,
                       xstr = xstr,
                       ystr = ystr,
                       Plot = False,
                       **kwargs,
                #      xerr=yerr0,
                      )
            except:
#                    traceback.print_exc()
                pass


    return ax



def PlotPaperVersionScatter(pathh_indx, 
                            xindx=0,
                            yindx=0,
                            **kwargs):
    verbose = True if 'verbose' in kwargs and kwargs['verbose'] == True else False
    #--- setup symbols
    colors = ['black','red','green','blue','cyan','brown','grey','magenta','orange','yellow']
    fillstyles=['white',None,'white',None,'white',None,'white',None,'white',None,'white',None,'white',None,'white',None]
    markers=['o','s','D','^','<','>','v']

    #--- plot
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    #
    ax.tick_params(labelsize=20,which='both',axis='both', top=True, right=True)
    #
    #
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'],fontsize=20 if not 'fontsize' in kwargs else kwargs['fontsize'])
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'],fontsize=20 if not 'fontsize' in kwargs else kwargs['fontsize'])
    ax.tick_params(labelsize=20 if not 'fontsize' in kwargs else kwargs['fontsize'])
    #
    Xdata=[]
    Ydata=[]
    for mg, marker, color, fillstyle in list(zip( [ 
                                        'FeNi',
                                           'CoNiFe',
                                           'CoNiCrFe',
                                            'CoCrFeMn',
                                            'CoNiCrFeMn',
                                              'Co5Cr2Fe40Mn27Ni26'
                                        ],markers, colors, fillstyles )):    
#         try:
#             os.system('rm crMuRho_concat.txt')
#         except:
#             pass
#         try:
#             os.system('rm hmin.txt')
#         except:
#             pass
        #--- extract data to be plotted
#        print('modify path!')
        
        #--- loop over time
        for itimee in [0]:
            
            
            #--- loop over realizations
            for irun in [0,1,2]: #range(3):
                pathh = { 
                          0:'%s/PairCrltnT300/%s/itime%s/Run%s'%(os.getcwd(),mg,itimee,irun),
                          1:'%s/VorAnlT300/%s/Run%s'%(os.getcwd(),mg,irun),
                          2:'%s/D2minAnalysisT300/%s/Run%s'%(os.getcwd(),mg,irun),
                        }[pathh_indx]
                file_dict = {
                          0:'%s/gr.txt'%pathh,
                          1:'%s/icoPercentageWithStrain.txt'%pathh,
                          2:'%s/d2min_gamma.txt'%pathh,
                          3:'%s/rc_d2min.txt'%pathh,
                          4:'%s/crsD2minRhoSro.txt'%pathh,
                          5:'%s/NegativeModulus.txt'%pathh,
                          6:'%s/YieldDelta.txt'%pathh,
                          7:'%s/alpha.txt'%pathh,
                        }
                file0 = file_dict[xindx] #--- xdata
                file1 = file_dict[yindx] #--- ydata
                xdata = [np.nan]
                ydata = [np.nan]
                if os.path.isfile( file0 ) and os.path.isfile( file1 ):
                    if verbose:
                        print('file0=',file0)
                        print('file1=',file1)
#                    os.system('cat %s >> crMuRho_concat.txt' %file0 )
#                    os.system('cat %s >> mu_size.txt'%file1)
                
                    #--- load data
                    sarr = np.loadtxt(file0)
                    sarr2nd = np.loadtxt(file1)
#                     pdb.set_trace()
                
                    try:
                        xdata = sarr[:,kwargs['colx']]
                    except: #--- only one row!
                        xdata = np.array(sarr[kwargs['colx']])
                    try:
                        ydata = sarr2nd[:,kwargs['coly']]
                    except:
                        ydata = np.array(sarr2nd[kwargs['coly']])
                        
                    #--- process
                    if xindx == 5 or xindx == 7: #--- hmin
                        xdata *= -1 
                    if yindx == 5 or yindx == 5: 
                        ydata *= -1
                    #    
                    #
                    #
                    if xindx == 4: #--- cxy
                        #--- fliter
                        #ydata = np.mean(ydata[-1-8:-1]) #--- mean of the last 10 points
                        #--- %95 quantile 
                        xdata.sort()
                        n=len(xdata)
                        xdata = xdata[int(.05*n):][0]
                    if yindx == 4: #--- cxy
                        #--- fliter
                        #ydata = np.mean(ydata[-1-8:-1]) #--- mean of the last 10 points
                        #--- %95 quantile 
                        ydata.sort()
                        n=len(ydata)
                        ydata = ydata[int(.05*n):][0]
                    #
                    #
                    #
                    if xindx == 3: #--- w_sb
                        xdata = sarr[:,1:3] #--- 2nd and 3rd cols are widths
                        indices = np.any(~np.isnan(xdata),axis=1)
                        xdata = xdata[~np.isnan(xdata)]
                        xdata = np.array([xdata[-1]]) #--- lastpoint
                    if yindx == 3: #--- w_sb
                        ydata = sarr2nd[:,1:3] #--- 2nd and 3rd cols are widths
                        indices = np.any(~np.isnan(ydata),axis=1)
#                        pdb.set_trace()
                        ydata = ydata[~np.isnan(ydata)]
                        ydata = np.array([ydata[-1]])

                    if verbose:
                        print('xdata=',xdata,'ydata=',ydata)
                #--- append data
                if type(xdata) != type([]): #--- avoid nans
                    Xdata.append(xdata)
                    Ydata.append(ydata)

                #--- plot different realizations
                #--- graph-related attributes
                attrs={ 'color':color,
                        'markersize':10 if not 'markersize' in kwargs else kwargs['markersize'] ,
                        'marker':marker,
                        'markerfacecolor':fillstyle,
                        'markeredgecolor':'black' if not fillstyle else None,
                        'label':'%s/irun%s'%(mg,irun),
                       'fmt':'.',
                      }
                # plot ydata
                if verbose:
                    print('kwargs=',kwargs)
                PltErr(xdata,ydata, 
                       yerr=None, #0*ydata,      
                       ax = ax,
                       attrs = attrs,
                       Plot = False,
                       **kwargs,
#                           xerr=yerr0,
                      )
            
    #
#     LOGY = True if ('yscale' in kwargs and kwargs['yscale'] == 'log') else False
#     LOGX = True if ('xscale' in kwargs and kwargs['xscale'] == 'log') else False
#     PutMinorTicks(ax, LOGX=LOGX,LOGY=LOGY)
#     #
#     DrawFrame(ax, 0.2,0.09,0.15,0.06,
#               kwargs['borderwidth'] if 'borderwidth' in kwargs else 0.01,
#               LOG_Y=LOGY,LOG_X=LOGX)
#     #
#     if 'legend' in kwargs and kwargs['legend']:
# #        fig.legend(bbox_to_anchor=(1.5,.8))
#         ax.legend(frameon=False, fontsize=12, markerscale=0.8,handlelength=1.3,handletextpad=0.4)
#     #
#     plt.savefig(kwargs['title'] if 'title' in kwargs else 'junk.png',dpi=2*75,bbox_inches='tight',
#                 pad_inches=0.0)
#     plt.show()
    return ax

def GetMetrics(pathh_indx,file0_indx,**kwargs):
    #--- setup symbols
    colors = ['black','red','green','blue','cyan','brown','grey','magenta','orange','yellow']
    fillstyles=['white',None,'white',None,'white',None,'white',None,'white',None,'white',None,'white',None,'white',None]
    markers=['o','s','D','^','<','>','v']

    #--- plot
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    #
    ax.tick_params(labelsize=20,which='both',axis='both', top=True, right=True)
    #
    #
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'],fontsize=20)
    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'],fontsize=20)
    ax.tick_params(labelsize=20)
    #
    for mg, marker, color, fillstyle in list(zip( [ 
                                        'FeNi',
                                          'CoNiFe',
                                           'CoNiCrFe',
                                            'CoCrFeMn',
                                            'CoNiCrFeMn',
                                              'Co5Cr2Fe40Mn27Ni26'
                                        ],markers, colors, fillstyles )):    

        #--- extract data to be plotted
#        print('modify path!')
        
        #--- loop over time
        for itimee in [0]:
            
            #--- store ensemble average
            Xdata = None
            Ydata = None
            
            #--- loop over realizations
            for irun in [0,1,2]: #range(3):
                pathh = { 
                          0:'%s/PairCrltnT300/%s/Run%s'%(os.getcwd(),mg,irun),
                          1:'%s/VorAnlT300/%s/Run%s'%(os.getcwd(),mg,irun),
                          2:'%s/D2minAnalysisT300/%s/Run%s'%(os.getcwd(),mg,irun),
                        }[pathh_indx]
                file0 = {
                          0:'%s/gr.txt'%pathh,
                          1:'%s/icoPercentageWithStrain.txt'%pathh,
                          2:'%s/d2min_gamma.txt'%pathh,
                          3:'%s/rc_d2min.txt'%pathh,
                          4:'%s/crsD2minRhoSro.txt'%pathh,
                          5:'%s/NegativeModulus.txt'%pathh, 
                          6:'%s/YieldDelta.txt'%pathh, 
                          7:'%s/gr_peak_gamma.txt'%pathh, 
                        }[file0_indx]
                
                xdata = [np.nan]
                ydata = [np.nan]
                if os.path.isfile( file0 ): # and os.path.isfile( file1 ):
#                    os.system('cat %s >> crMuRho_concat.txt' %file0 )
#                    os.system('cat %s >> mu_size.txt'%file1)
                
                    #--- load data
                    sarr = np.loadtxt(file0)
#                    sarr2nd = np.loadtxt(file1)
                    try:
                    #--- assign x, y
#                    print(sarr.shape)
                        if file0_indx == 5:
                            xdata = np.array([sarr[0]])
                            ydata = np.array([-sarr[1]])
                        elif file0_indx == 6:
                            xdata = np.array([sarr[0]])
                            ydata = np.array([sarr[1]])
                        else:
                            xdata = sarr[:,0] #--- x
                            ydata = sarr[:,1] #--- y

                        #--- only for shear band width
                        if file0_indx == 3:
                            ydata2 = sarr[:,2] #--- y1
                            #--- remove nan's
                            indices = np.any(~np.isnan(np.c_[ydata,ydata2]),axis=1) #--- true:either x or y is not nan
                            ydata = np.array(list(map(lambda x:np.mean(x[~np.isnan(x)]),np.c_[ydata,ydata2][indices]))) #--- take mean                    
                            xdata = xdata[indices]
                            assert xdata.shape == ydata.shape
                    
                    except:
                        if verbose: print('error while reading %s'%file0)
                        pass

#                    ydata /= xdata**2 #!!!!!!!!!!!!!comment 
#                    yerr = sarr[:,2]
#                    pdb.set_trace()

                    #--- ensemble average
                    try:
                        Xdata = np.concatenate((Xdata,xdata),axis=0) 
                        Ydata = np.concatenate((Ydata,ydata),axis=0) 
                    except:
                        Xdata = xdata.copy()
                        Ydata = ydata.copy()
#                    print(Xdata.shape)

                    #--- plot different realizations
                if not 'PlotMean' in kwargs or not kwargs['PlotMean']:
                    #--- graph-related attributes
                    attrs={ 
                            'color':color,
                            'markersize':10,
                            'marker':marker,
                            'markerfacecolor':fillstyle,
                            'markeredgecolor':'black' if not fillstyle else None,
                            'label':'%s/irun%s'%(mg,irun),
                           'markevery':max(1,int(len(xdata)/10.0)),
                          }
                    # plot ydata
                    PltErr(xdata,ydata, 
                           yerr=None, #0*ydata,      
                           ax = ax,
                           attrs = attrs,
                           Plot = False,
                           **kwargs,
#                           xerr=yerr0,
                          )
                    #--- plot ydata2
#                         PltErr(xdata,ydata2, 
#                                yerr=None, #0*ydata,      
#                                ax = ax,
#                                attrs = attrs,
#                                Plot = False,
#                                **kwargs,
#     #                           xerr=yerr0,
#                               )

        #
        if 'PlotMean' in kwargs and kwargs['PlotMean']:
            #--- take average
            Xbin, Ybin, Yerr = BinData(Xdata,Ydata)
            attrs={ 'color':color,
                    'markersize':10,
                    'marker':marker,
                    'markerfacecolor':fillstyle,
                    'markeredgecolor':'black' if not fillstyle else None,
                    'label':'%s'%(mg),
                   'markevery':int(len(Xbin)/10.0),
                   'errorevery':int(len(Xbin)/10.0),
                   'markeredgewidth':0.7,
                    'linewidth':1, 
                   'barsabove':None,
                   'capsize':5,
                   'capthick':1,
                   'elinewidth':1
                  }

            PltErr(Xbin,Ybin, 
                   yerr=Yerr,      
                   ax = ax,
                   attrs = attrs,
                   Plot = False,
                   **kwargs,
            #      xerr=yerr0,
                  )

        
    #
    LOGY = True if ('yscale' in kwargs and kwargs['yscale'] == 'log') else False
    LOGX = True if ('xscale' in kwargs and kwargs['xscale'] == 'log') else False
    PutMinorTicks(ax, LOGX=LOGX,LOGY=LOGY)
    #
    DrawFrame(ax, 0.2,0.09,0.15,0.06,
              kwargs['borderwidth'] if 'borderwidth' in kwargs else 0.01,
              LOG_Y=LOGY,LOG_X=LOGX)
    #
    if 'legend' in kwargs and kwargs['legend']:
#        fig.legend(bbox_to_anchor=(1.5,.8))
        ax.legend(frameon=False, fontsize=12, markerscale=0.8,handlelength=1.3,handletextpad=0.4)
    #
#    plt.savefig(kwargs['title'] if 'title' in kwargs else 'junk.png',dpi=2*75,bbox_inches='tight',
#                pad_inches=0.0)
    plt.show()
    return ax




def ReturnShapeFunc(x,df):
#    eps1=0.0;eps2=0.0
#    n = len(x)
#    x2 = x * x
#    x3 = x * x2
#    x4 = (x-eps1)**3
#    x5 = (x-eps2)**3
    transformed_x = patsy.bs(x, df=df,degree =3, include_intercept=False)
    return transformed_x 

def TrainModel(x,y):
    reg = LinearRegression().fit(x, y )
    return reg
def Validate(reg, x,y,deg_f=range(0,40,10),cv_samples=10):
    mse={}
    for df in deg_f:
#    for df in map(int,np.logspace(0,12,20,base=2,endpoint=True)):
        try:
            transfrm = ReturnShapeFunc(x,df)
            #--- vadidate    
            scores = cross_validate(reg,  transfrm, y, cv=cv_samples,
                                         scoring=('r2', 'neg_mean_squared_error'),
                                         return_train_score=True)
            mse[df] = np.mean(scores['train_neg_mean_squared_error'])
        except:
            continue
    return mse

def YieldStress(x,y_pred):
    #--- max. stress
    y_pred_arr = np.array(y_pred.flatten())
    indx = np.arange(len(y_pred_arr))[y_pred_arr==np.max(y_pred_arr)][0]
    ey =  x[indx]
#    print('ey=',ey)

    #--- yield stress
    sy =  y_pred_arr[indx]
#    print('sy=',sy)
    assert sy > 0.0

    #--- max. slope
    ydiff = np.gradient(y_pred_arr, x)
    indx = np.arange(len(x))[ydiff==np.min(ydiff[x>ey])][0]
    em = x[indx]
    sm =  y_pred_arr[indx]
    dsm = ydiff[x>=em][0]
#    print('em=',em)

    #--- flow stress
#    ydiff = np.gradient(y_pred_arr, x)
    indx = np.arange(len(x))[np.abs(ydiff)==np.min(np.abs(ydiff)[x>em])][0]
    ef = x[indx]
    sf =  y_pred_arr[indx]
    
    #--- modulus
    xdata = x[x<0.5*ey]
    ydata = y_pred[x<0.5*ey]
    xdata -= np.mean(xdata)
    ydata -= np.mean(ydata)
    gmodu = np.polyfit(xdata, ydata, 1)[0] #ydiff[0]

    return (ey,sy), (em,sm,dsm), (ef,sf), gmodu

def YieldStress2nd(x,y_pred):
    #--- max. stress
    y_pred_arr = np.array(y_pred.flatten())
    indx = np.arange(len(y_pred_arr))[y_pred_arr==np.max(y_pred_arr)][0]
    ey =  x[indx]
#    print('ey=',ey)

    #--- yield stress
    sy =  y_pred_arr[indx]
#    print('sy=',sy)
    assert sy > 0.0



    #--- flow stress
    ydiff = np.gradient(y_pred_arr, x)
    indx = np.arange(len(x))[np.abs(ydiff)==np.min(np.abs(ydiff)[x>ey])][0]
    ef = x[indx]
    sf =  y_pred_arr[indx]


    #--- max. slope
#    ydiff = np.gradient(y_pred_arr, x)
    indx = np.arange(len(x))[ydiff==np.min(ydiff[x<ef])][0]
    em = x[indx]
    sm =  y_pred_arr[indx]
    dsm = ydiff[x>=em][0]
#    print('em=',em)

    #--- modulus
    xdata = x[x<0.5*ey]
    ydata = y_pred[x<0.5*ey]
    xdata -= np.mean(xdata)
    ydata -= np.mean(ydata)
    gmodu = np.polyfit(xdata, ydata, 1)[0] #ydiff[0]
    
    
    return (ey,sy), (em,sm,dsm), (ef,sf), gmodu

def BinData(Xdata,Ydata,nbins=1024, scale='linear'):
        if scale == 'linear':
            bins =nbins
        elif scale=='log':
            filtr=Filtr(Xdata,0.0,np.inf)
            xmin = np.log10(Xdata[filtr].min())
            xmax = np.log10(Xdata[filtr].max())
            bins = np.logspace(xmin,xmax,nbins)
        ysum, edges = np.histogram(Xdata,weights=Ydata,bins=bins)
        ysum_sq, edges = np.histogram(Xdata,weights=Ydata*Ydata,bins=bins)
        xsum, edges = np.histogram(Xdata,weights=Xdata,bins=bins)
        xcount, edges = np.histogram(Xdata,bins=bins)
        
        xsum=xsum[xcount>0]
        ysum=ysum[xcount>0]
        ysum_sq=ysum_sq[xcount>0]
        xcount=xcount[xcount>0]
        #--- std. deviation
        ymean = ysum / xcount
        ysq_mean = ysum_sq / xcount
        ystd = ( ysq_mean - ymean * ymean ) ** 0.5
        
        return xsum / xcount, ysum / xcount, ystd / xcount ** 0.5


    


def PltBitmapWithScatter( value, xyScatter,
              xlim = (-0.5,0.5), ylim = (-0.5,0.5),
              **kwargs
             ):
        
    val = value.copy()
    #--- z-score
    if 'zscore' in kwargs and kwargs['zscore']:
        val -= np.mean(val)
        val /= np.std(val)
        val[val>2.0]=1.0
        val[val<-2.0]=-1.0
    #--- plot
    (mgrid,ngrid) = val.shape
    center = (ngrid/2,mgrid/2)
    #
#    aspect = (ylim[1]-ylim[0])/(xlim[1]-xlim[0])
    fig = plt.figure(figsize=(4,4)) #*aspect))
    ax = fig.add_subplot(111)
    #---
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    #
    fontsize = kwargs['fontsize'] if 'fontsize' in kwargs else 16
#    print('fontsize=',fontsize)
    #
    ax.set_xlabel(kwargs['xlabel'] if 'xlabel' in kwargs else '',fontsize=fontsize)
    ax.set_ylabel(kwargs['ylabel'] if 'ylabel' in kwargs else '',fontsize=fontsize)
    #
    ax.tick_params(labelsize=fontsize,which='both',axis='both', top=False, right=False)
    #
    if 'ticklabels' in kwargs:
        ax.axes.xaxis.set_visible(kwargs['ticklabels'])
        ax.axes.yaxis.set_visible(kwargs['ticklabels'])
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(kwargs['ticklabels'])
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(kwargs['ticklabels'])
    #
    pos = ax.imshow(val.real,cmap='bwr',
                     extent=(xlim[0],xlim[1],ylim[0],ylim[1]),origin='lower')# ,vmin=-.6, vmax=.6)
    
    ax.scatter(xyScatter[:,0],xyScatter[:,1],
           alpha=1,color=kwargs['color'] if 'color' in kwargs else 'yellow',
               marker='o',
               s=kwargs['s'] if 's' in kwargs else 4)
    #
    PutMinorTicks(ax)
    #
    if 'colorbar' in kwargs and kwargs['colorbar']:
        fig.colorbar( pos, fraction = 0.04)
    if 'DrawFrame' in kwargs: 
        DrawFrame(ax, *kwargs['DrawFrame'])
    if 'title' in kwargs:
        plt.savefig(kwargs['title'],dpi=kwargs['dpi'] if 'dpi' in kwargs else 75,bbox_inches='tight',pad_inches=0.0)
    plt.show()
    


def PlotDminVor( lmpData, lmpDmin, times,
                   title='d2min',
               ):
    
    for itime in sorted(times):
#         pdb.set_trace()
        #--- dmin
        (xlin, ylin, zlin), junk, d2intrp = Intrp(lmpDmin, 
                                            times = [itime],
                                            time0 =0,
                                            Plot = None,
                                           )
#         pdb.set_trace()
        #--- reshape value
        nx,ny,nz = len(xlin), len(ylin),len(zlin) 
        value = np.c_[d2intrp[itime].d2min].reshape(((ny,nx,nz)))



        #--- icosahedra
        vor = lp.Atoms( **lmpData.coord_atoms_broken[itime].to_dict(orient='list') )
        
        #--- map
        box = lp.Box( BoxBounds = lmpData.BoxBounds[itime] )
        box.BasisVectors( AddMissing = np.array([0.0,0.0,0.0] ))
        #
        mapp = lp.Map( vor, box ) 
        mapp.ChangeBasis()
        mapp.Set( vor ) #--- atoms: add mapped xyz

        #--- filter
        indices = GetFullIcosahedra( vor )
        dz = zlin[-1]-zlin[-2]
        dz *= 1.25
        lz = zlin[-1]-zlin[0]
        zc=0.5*(zlin[0]+zlin[-1])+10*dz
        zlo=(zc-dz-zlin[0]) % lz + zlin[0]
        zhi=(zc+dz-zlin[0]) % lz + zlin[0]
        
#        print(zlo,zc,zhi)
#        pdb.set_trace()
        indxx= np.all([vor.zm > zlo, 
                       vor.zm < zhi
                      ],
                        axis=0)
        indxx = np.all([indxx,indices],axis=0)

                #--- xy plane
        
        val = value[:,:,int(nz*(zc-zlo)/lz)].copy()

        #
        PltBitmapWithScatter(val, 
                      np.c_[np.array(vor.xm)[indxx],np.array(vor.ym)[indxx]],
                      xlabel = '', ylabel = '',
    #                  xlim=VectorNorm[0]*np.array([0.0,1.0]),ylim=VectorNorm[1]*np.array([0.0,1.0]),
                      xlim=np.array([xlin[0],xlin[-1]]),ylim=np.array([ylin[0],ylin[-1]]),
                      zscore = True,
                      frac = 1.0, #--- plot a patch
                      title = '%s%s.png'%(title,itime),
 #                     colorbar=None,
 #                     ticklabels = None,
                      color='black',
 #                     DrawFrame=[0.2,0.09,0.15,0.06,0.04],
                      )

def PltErr( xdata, ydata, 
            yerr = None,
            xstr = '',
            ystr = '',
            Plot = True,
            **kwargs,
            ):
    fontsize=kwargs['fontsize'] if 'fontsize' in kwargs else 20
    if not 'ax' in kwargs:
        fig = plt.figure( figsize = (4,4))
        ax = fig.add_subplot(111)
    else:
        ax = kwargs['ax']
        if 'twinx' in kwargs and kwargs['twinx']:
                ax = kwargs['ax'].twinx()
    #--- setting   
    ax.set_xlabel(xstr,fontsize=fontsize)
    ax.set_ylabel(ystr,fontsize=fontsize)
    ax.tick_params(labelsize=fontsize,which='both',axis='both', top=True, right=True)
    #
    xerr = kwargs['xerr'] if 'xerr' in kwargs else None 
#
    if 'attrs' in kwargs:
        ax.errorbar( xdata, ydata,yerr = yerr, xerr = xerr, **kwargs['attrs'])
        if 'fill_between' in kwargs and kwargs['fill_between']:   
            ax.fill_between(xdata, ydata-yerr, ydata+yerr)

    else:
        ax.errorbar( xdata, ydata,yerr = yerr, xerr = xerr, fmt='-o',label=r'$x$')       
    #--- plot
    #
#    ax.plot(ax.axis()[:2],[0.0,0.0],'-.',color='black')
    #
    if 'ylim' in kwargs:
        ylim = kwargs['ylim'] 
        ax.set_ylim(ylim)
    if 'xlim' in kwargs:
        xlim = kwargs['xlim'] 
        ax.set_xlim(xlim)
    #
    if 'xscale' in kwargs: 
        ax.set_xscale(kwargs['xscale'])
    if 'yscale' in kwargs: 
        ax.set_yscale(kwargs['yscale'])
    #
    if 'xticks' in kwargs:
        ax.set_xticks(list(map(float,kwargs['xticks'][1])))
        ax.set_xticklabels(list(map(lambda x:'$%s$'%x,kwargs['xticks'][0])))
    #
    if 'yticks' in kwargs:
        ax.set_yticks(list(map(float,kwargs['yticks'][1])))
        ax.set_yticklabels(list(map(lambda x:'$%s$'%x,kwargs['yticks'][0])))
        
    #
    LOGY = True if ('yscale' in kwargs and kwargs['yscale'] == 'log') else False
    LOGX = True if ('xscale' in kwargs and kwargs['xscale'] == 'log') else False
    PutMinorTicks(ax, LOGX=LOGX,LOGY=LOGY)
    #
    if 'DrawFrame' in kwargs: 
        DrawFrame(ax, *kwargs['DrawFrame'],LOG_Y=LOGY,LOG_X=LOGX)
    #
    if 'legend' in kwargs:# and :
        plt.legend(**kwargs['legend'])#frameon=False,fontsize=fontsize)
    if 'title' in kwargs: #Plot:
        plt.savefig(kwargs['title'],dpi=300 if not 'dpi' in kwargs else kwargs['dpi'],bbox_inches='tight', 
                    pad_inches=0.0)
    if Plot:
        plt.show()
    #
    
    
#    if not 'ax' in kwargs:
    return ax
       

def PltScatter( xdata, ydata, 
            yerr = None,
            xstr = '',
            ystr = '',
            Plot = True,
            **kwargs,
            ):
    fontsize=kwargs['fontsize'] if 'fontsize' in kwargs else 20
    if not 'ax' in kwargs:
        fig = plt.figure( figsize = (4,4))
        ax = fig.add_subplot(111)
    else:
        ax = kwargs['ax']
        if 'twinx' in kwargs and kwargs['twinx']:
                ax = kwargs['ax'].twinx()
    #--- setting   
    ax.set_xlabel(xstr,fontsize=fontsize)
    ax.set_ylabel(ystr,fontsize=fontsize)
    ax.tick_params(labelsize=fontsize,which='both',axis='both', top=True, right=True)
    #
    xerr = kwargs['xerr'] if 'xerr' in kwargs else None 
#
    if 'attrs' in kwargs:
        ax.scatter( xdata, ydata, **kwargs['attrs'])
        if 'fill_between' in kwargs and kwargs['fill_between']:   
            ax.fill_between(xdata, ydata-yerr, ydata+yerr)

    else:
        ax.errorbar( xdata, ydata,yerr = yerr, xerr = xerr, fmt='-o',label=r'$x$')       
    #--- plot
    #
#    ax.plot(ax.axis()[:2],[0.0,0.0],'-.',color='black')
    #
    if 'ylim' in kwargs:
        ylim = kwargs['ylim'] 
        ax.set_ylim(ylim)
    if 'xlim' in kwargs:
        xlim = kwargs['xlim'] 
        ax.set_xlim(xlim)
    #
    if 'xscale' in kwargs: 
        ax.set_xscale(kwargs['xscale'])
    if 'yscale' in kwargs: 
        ax.set_yscale(kwargs['yscale'])
    #
    if 'xticks' in kwargs:
        ax.set_xticks(list(map(float,kwargs['xticks'][1])))
        ax.set_xticklabels(list(map(lambda x:'$%s$'%x,kwargs['xticks'][0])))
    #
    if 'yticks' in kwargs:
        ax.set_yticks(list(map(float,kwargs['yticks'][1])))
        ax.set_yticklabels(list(map(lambda x:'$%s$'%x,kwargs['yticks'][0])))
        
    #
    LOGY = True if ('yscale' in kwargs and kwargs['yscale'] == 'log') else False
    LOGX = True if ('xscale' in kwargs and kwargs['xscale'] == 'log') else False
    PutMinorTicks(ax, LOGX=LOGX,LOGY=LOGY)
    #
    if 'DrawFrame' in kwargs: 
        DrawFrame(ax, *kwargs['DrawFrame'],LOG_Y=LOGY,LOG_X=LOGX)
    #
    if 'legend' in kwargs and kwargs['legend']:
        plt.legend(frameon=False,fontsize=fontsize)
    if 'title' in kwargs: #Plot:
        plt.savefig(kwargs['title'],dpi=300 if not 'dpi' in kwargs else kwargs['dpi'],bbox_inches='tight', 
                    pad_inches=0.0)
    if Plot:
        plt.show()
    #
    
    
#    if not 'ax' in kwargs:
    return ax
        
def FilterDataFrame(df,key='id',val=[1,2,3]): #,out='C66'):
    tmp0 = df.set_index(key,drop=True,append=False).loc[val] 
    return tmp0.reset_index() #.reindex(range(len(tmp0)))

def CrssCrltn(x,y):        
    x-=np.mean(x)
    x/=np.std(x)

    y-=np.mean(y)
    y/=np.std(y)

    return (np.sum(x*y)/len(x))


def ScatterXY( vor, d2min,
                  Plot = True, PLOT_AVERAGE = True, title='scatterD2minRho.png',
                 GetAx = False,
                nbins_per_decade = 8,
                 **kwargs,
               
                ):

            #
    xscale = 'log' if not 'xscale' in kwargs else kwargs['xscale']
    yscale = 'log' if not 'yscale' in kwargs else kwargs['yscale']

    #--- cross crltn.
    x = np.array(vor.tmp)
    y = np.array(d2min.tmp)
    #---
    if yscale == 'log':
        y=np.log10(y) #y[x>0]
    if xscale == 'log':
        x=np.log10(x) #y[x>0]
#         x=x[x>0]
        
    if 'zscore' in kwargs and kwargs['zscore']:
        x=Zscore(x)
        y=Zscore(y)

    crs = CrssCrltn( x,y )
        
    if Plot:
        #--- plot
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        #
        ax.scatter(x,y,marker='x',alpha=0.08)
        #
        if PLOT_AVERAGE:
            xb, yb, err_x, err_y = GetBinnedAverage( x, y, nbins_per_decade = nbins_per_decade )
            nevery = int(np.ceil(len(xb)/16.0))
#            pdb.set_trace()
#            lower_errory =  yb-err_y #10**yb*(1-10**(-err_y))
#            upper_errory =  yb+err_y #10**yb*(10**(+err_y)-1)
#            asymmetric_errory = np.array(list(zip(lower_errory, upper_errory))).T
            #
#            err_x = np.zeros(len(xb)) if 'xerr' in kwargs and not kwargs['xerr'] else err_x
#            print(err_x)
#            lower_errorx =  xb-err_x  #10**xb*(1-10**(-err_x))
#            upper_errorx =  xb+err_x #10**xb*(10**(+err_x)-1)
#            asymmetric_errorx = np.array(list(zip(lower_errorx, upper_errorx))).T
            
#             ax.errorbar(x=xb,y=yb,yerr=asymmetric_errory,xerr=asymmetric_errorx,
#                 fmt='-o', markersize=8,color='red',markerfacecolor='white',markeredgewidth=0.7,
#                         markeredgecolor='red', markevery=nevery, errorevery=nevery,
#                 linewidth=.5, barsabove=None,capsize=5,capthick=1,elinewidth=1)
            ax.errorbar(x=xb,y=yb,yerr=err_y,xerr=err_x,
                fmt='-o', markersize=8,color='red',markerfacecolor='white',markeredgewidth=0.7,
                        markeredgecolor='red', markevery=nevery, errorevery=nevery,
                linewidth=.5, barsabove=None,capsize=5,capthick=1,elinewidth=1)

#                 def func(x,a,b):
#                     return a*x+b
#                 z = np.polyfit(xb, yb, 1)

#                 label = r'$%2.1e.x+%2.1e$'%(z[0],z[1])
            
        if 'xlim' in kwargs:
            ax.axis(kwargs['xlim'])
#        makeTickLabels(ax, x, y)
        #
#        ax.set_yscale(yscale)
#        ax.set_xscale(xscale)
        #
        if 'axisLabels' in kwargs:
            (xlabel, ylabel) = kwargs['axisLabels']
            ax.set_ylabel(ylabel,fontsize=16)
            ax.set_xlabel(xlabel,fontsize=16)
            #
#            ax.set_title(r'cxy=%3.2e'%crs,fontsize=16)
            #
#        ax.legend(frameon=False, fontsize=12)
        #
        ax.tick_params(labelsize=20,which='both',axis='both', top=True, right=True)
        
        DrawFrame(ax, 0.2,0.09,0.15,0.06,0.04,
                  LOG_X=True if xscale == 'log' else False,
                  LOG_Y=True if yscale == 'log' else False) 
        #
        plt.savefig(title,dpi=2*75,bbox_inches='tight',pad_inches=0.0)
        plt.show()
    
    if GetAx:
        return crs, ax
    else:
        return crs


def GetFullIcosahedra( vor ):
    dff=pd.DataFrame(vor.__dict__)
    indices = np.all([  dff['VoronoiIndex1']==0.0,dff['VoronoiIndex2']==0.0,dff['VoronoiIndex3']==0.0,\
                    dff['VoronoiIndex4']==0.0,dff['VoronoiIndex5']==12.0,dff['VoronoiIndex6']==0.0,\
                    dff['VoronoiIndex7']==0.0,dff['VoronoiIndex8']==0.0,dff['VoronoiIndex9']==0.0],axis=0)
#     indices = np.all([  dff['VoronoiIndex3']==0.0,dff['VoronoiIndex4']==0.0,dff['VoronoiIndex5']==12.0,\
#                     dff['VoronoiIndex6']==0.0],axis=0)
#    pdb.set_trace()
#    vor_df = pd.DataFrame(vor.__dict__)
        
#    df = pd.DataFrame(np.concatenate((np.c_[vor_df], np.c_[indices]), axis=1, dtype=np.object),
#                      columns=list(vor_df.keys())+['IcoTrue'])

    return indices #lp.Atoms(**df.to_dict(orient='series'))


def Get2dSlice( value, zlin, zc, nzll=[0] ):
        #--- get xy plane
#    zc=0.5*(zlin[0]+zlin[-1])
    dz = zlin[-1]-zlin[-2]
    lz = zlin[-1]-zlin[0]
    nz = len(zlin)
    #
    zz = zc #zlin[-1] #zc #zlin[-1] #--- arbitrary plane
    nzz=int(nz*(zz-zlin[0])/lz)
#    print(nzz,nz)
    if nzz == nz: nzz -= 1
    val = value[:,:,nzz].copy()
    #
    nzll[0] = nzz
    return val

def Get2dSliceScattered( xm, ym, zm, 
                        zlin, zc ):    
    dz = zlin[-1]-zlin[-2]
    lz = zlin[-1]-zlin[0]
    nz = len(zlin)
    dz *= 2.0
    zz=zc
    #
    zlo=(zz-dz-zlin[0]) % lz + zlin[0]
    zhi=(zz+dz-zlin[0]) % lz + zlin[0]
    indxx= np.all([zm > zlo, 
                   zm < zhi
                  ],
                    axis=0)
    if zhi < zlo: #--- periodic boundaries
        indxx= np.any([zm > zlo, 
           zm < zhi
          ],
            axis=0)
    return np.c_[xm,ym][indxx]        
    
def GetInterpolatedData( x,y,z, 
                         xv, yv, zv,
                       ):
    r2nd=4.8 #--- g(r)
    dx = xv[1]-xv[0]
    dy = yv[1]-yv[0]
    dz = zv[1]-zv[0]
    dr = (dx*dx+dy*dy+dz*dz)**0.5
    sigma = int(r2nd/dr)
    heatmap, edges = np.histogramdd( np.c_[y, x, z],
                                    bins=[np.append(yv,yv[-1]+dy),
                                          np.append(xv,xv[-1]+dx),
                                          np.append(zv,zv[-1]+dz)],
                                    normed=True)

    print('(nx,ny,nz)=', heatmap.shape )
    
    heatmap *= len( x )
    heatmap = gaussian_filter( heatmap, sigma = sigma )
    return heatmap



def DensityMap(atoms0, filtr, xv, yv,zv):    
    #--- scattered points
    xs = np.array(atoms0.xm)[filtr] #--- must be reference frame!
    ys = np.array(atoms0.ym)[filtr]
    zs = np.array(atoms0.zm)[filtr]

    (nx,ny,nz)=list(map(len,[xv,yv,zv]))
    #--- density
    heatmap=GetInterpolatedData(xs,ys,zs,
                                xv, yv, zv,
                               )
    #--- 2d slice
    heatmap2d = Get2dSlice( heatmap, zv, 
                            zv[-1] )
    #--- plot scattered data
    xy = Get2dSliceScattered( xs, ys, zs, 
                        zv, zv[-1] )
    #--- bitmap
    PltBitmapWithScatter(
          heatmap2d, #transpose? 
          xy,
          xlabel = r'$x$(\r{A})', ylabel = r'$y$(\r{A})',
          xlim=np.array([xv.min(),xv.max()]),ylim=np.array([yv.min(),yv.max()]),
          zscore = True,
          frac = 1.0, #--- plot a patch
          title = 'kdeRho.png',
#          colorbar=True,
#          ticklabels = True,
          color='black',
#                     DrawFrame=[0.2,0.09,0.15,0.06,0.04],
          )

    
    return heatmap

def Zscore(val):
    x=val.copy()
    x -= np.mean(x)
    x /= np.std(x)
    return x


def SetTickLabels(ax, **kwargs):
    fmt='%3.1f'
    if 'xaxis' in kwargs:
        tickLabels = kwargs['xaxis']
        ax.xaxis.set_ticklabels(['$%s$'%i for i in tickLabels])
        ax.xaxis.set_ticks(tickLabels)
    if 'yaxis' in kwargs:
        tickLabels = kwargs['yaxis']
        ax.yaxis.set_ticklabels(['$%s$'%i for i in tickLabels])
        ax.yaxis.set_ticks(tickLabels)
        

 
        
def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
    
def gaussian_mixture( values, 
#                     times,
                     LABELS = False, 
                     PLOT = True):
#    thresh = {}
#    ratio = {}
#    pdb.set_trace()
    gsFitTotal = 1.0
    if 1:
#    for itime in sorted( times ):
    #--- convert data frema to array
        nij = values#[itime]
        X=np.log10(nij)
        X=np.array(X).reshape(len(X),1)

        try:
            gm_obj = skm.BayesianGaussianMixture(n_components=2, tol=1e-8, max_iter=10000,
                                                verbose=0, 
    #                                            random_state=0,
                                                init_params='random',
                                               )
            gm_obj.fit(X)

#            assert gm_obj.converged_

            #--- min(\mu0, \mu1) corresponds to trig. mode
            mean0 = gm_obj.means_[0][0]
            mean1 = gm_obj.means_[1][0]
            d     = { 'trigrd':min([mean0,0],[mean1,1])[1], 
                      'backgrd':max([mean0,0],[mean1,1])[1]}
            mean  = { 'trigrd' : gm_obj.means_[d['trigrd']][0],
                      'backgrd': gm_obj.means_[d['backgrd']][0] }
            sigma = { 'trigrd' : gm_obj.covariances_[d['trigrd']][0]**0.5,
                      'backgrd': gm_obj.covariances_[d['backgrd']][0]**0.5 }
            delta = { 'trigrd' : gm_obj.weights_[d['trigrd']],
                      'backgrd': gm_obj.weights_[d['backgrd']] }
    #        print(sigma['backgrd']/mean['backgrd'])
            #--- plot scatter
            nij_red = nij[gm_obj.predict(X)==0]
            nij_blue = nij[gm_obj.predict(X)==1]
            #ratio[itime] =  sigma['backgrd']/mean['backgrd'] #
            ratio = 1.0 * len(nij_blue) / (len(nij_red)+len(nij_blue))
    #        print(len(nij_red)+len(nij_blue))
    #        prob = np.random.rand(len(X))
    #         list_of_red = prob  <   gm_obj.predict_proba(X)[:, d['trigrd']] #--- assign modes based on the responsibilities \gamma
    #         list_of_blue = prob >= gm_obj.predict_proba(X)[:, d['trigrd']]
    #         nij_red = nij[list_of_red]
    #         nij_blue = nij[list_of_blue]

            #--- plot distributions
            edge_act, hist_act = DistNij(nij,normed=None, nbins_per_decade = 32)

            if PLOT:
                fig = plt.figure(figsize=(4,4))
                ax = fig.add_subplot(111)
    #            ax.set_yscale('log')
    #            ax.set_ylim(0.9,len(nij_red))#1e3)
                #
                ax.plot(edge_act,hist_act,'o',color='black',label='Total')
                #
                xv = edge_act
                ax.plot( xv, 
                        len(X)*(delta['trigrd']*gaussian(xv, mean['trigrd'], sigma['trigrd'])+
                                delta['backgrd']*gaussian(xv, mean['backgrd'], sigma['backgrd']))*(xv[1]-xv[0]), 
                        color='black')
                ax.plot( xv, len(X)*delta['trigrd']*gaussian(xv, mean['trigrd'], sigma['trigrd'])*(xv[1]-xv[0]),color='C0')
                ax.plot( xv, len(X)*delta['backgrd']*gaussian(xv, mean['backgrd'], sigma['backgrd'])*(xv[1]-xv[0]),color='red')
                if LABELS:
                    ax.set_xlabel(r'$log n_{ij}$')
                    ax.set_ylabel(r'$P(log n_{ij})$')
                ax.set_xlim(np.floor(np.min(edge_act)),np.ceil(np.max(edge_act)))
                #
                gsFitTotal = np.c_[ax.get_lines()[1].get_xdata(),ax.get_lines()[1].get_ydata()] #--- return data
            #


            #--- find decision boundary
            mu0,mu1,sigma0,sigma1,delta0,delta1=mean['trigrd'],mean['backgrd'],sigma['trigrd'],sigma['backgrd'],delta['trigrd'],delta['backgrd']
            def f(x): #,:
                return delta0*np.exp(-0.5*((x-mu0)/sigma0)**2)/sigma0 -\
                        delta1*np.exp(-0.5*((x-mu1)/sigma1)**2)/sigma1
            n_th = 0
            try:
                n_th = optimize.bisect(f, mu0-1*sigma0, mu1+1*sigma1)
                if PLOT:
                    ax.plot([n_th,n_th],ax.axis()[2:],'-.r') #--- vertical line
                thresh=10**n_th
            except:
                traceback.print_exc()
                pass
            if PLOT:
                fig.savefig('distnijGM.png',bbox_inches='tight',dpi=2*75)
                plt.show()

        except:
            traceback.print_exc()
            pass
    #     print gm_obj.predict_proba([[n_th]])
    return thresh, ratio, gsFitTotal

def DistNij(nij,normed=True, nbins_per_decade = 4, **kwargs):
#--- histogram
    nn=np.log10(nij)
    nmin=kwargs['nmin'] if 'nmin' in kwargs else nn.min()
    nmax=kwargs['nmax'] if 'nmax' in kwargs else nn.max()
    bins=np.linspace(nmin,nmax,int(nmax-nmin)*nbins_per_decade)
    hist, edge = np.histogram(nn,bins=bins,normed=normed)
    
    #--- accumulated histogram
    slist=np.array(nn)
    slist.sort()
    N = len(slist)
    d = histogramACCUMLTD( slist.tolist() )
    keys=list(d.keys())
    keys.sort()
    
    xx=[];yy=[]
    for ikey in keys:
        xx.append(d[ikey][0])
        yy.append(d[ikey][2])
    
    # ax2 = ax.twinx()
    # ax2.plot(xx,yy,
    #         linestyle='-', drawstyle='steps-post',color='red',
    #         linewidth=1.0) #--- accumulated
    # #ax2.set_xlim(-7,1)
    # #ax2.set_ylim(0,1200)
    
    # ax2.tick_params(axis='y',colors='red')
    # ax2.set_ylabel('$N(<n_{ij})$',color='red')
    
    return edge[:-1],hist

def valuesDict(d,keys):
    return list(map(d.get,sorted(keys)))


def histogramACCUMLTD( slist ):
    assert type( slist ) == type( [] ), 'arg must be a list. a %s is given!' %( type( slist ) )
    d = {}
    for item in slist:
        try:
            d[ item ] += 1
        except:
            d[ item ] = 1
    keys = list(d.keys())
    keys.sort()

    cdf = 0.0
    xi = min( slist ) - 1.0e-6
    xf = max( slist ) + 1.0e-6
    npoin = len( slist )
    adict = {}
    for ikey, index in zip( keys, range( sys.maxsize ) ):
        adict[ index ] = [ xi, ikey, cdf ]
        cdf += 1.0 * d[ ikey ] # / npoin
        xi = ikey
    adict[ index + 1 ] = [ xi, xf, cdf ]
    return adict

def GetStrain(lmpData, times, time0 ):
    ebulk = {}
    box0 = lp.Box( BoxBounds = lmpData.BoxBounds[time0], AddMissing = np.array([0.0,0.0,0.0] ))
    if 1:
#    for itime in times:
        itime = times
        box = lp.Box( BoxBounds = lmpData.BoxBounds[itime],AddMissing = np.array([0.0,0.0,0.0] ))
        #
        #--- bulk strain
        dx=box.CellVector[0,1]-box0.CellVector[0,1]
        l1=box.CellVector[1,1]
        ebulk = dx/l1
    return ebulk
        
    
def to_dict( df ):
    sdict = {}
    skeys = df.keys().to_list()
    for ikey in skeys:
        sdict[ikey] = df[ikey].to_list()

    return sdict

def isSane( AddColumns, columnList ):
    #--- avoid duplicates
    assert len( set( AddColumns ) ) == len( AddColumns ), 'duplicates in the list!'
    #--- assert column list is not already included
    n = len( AddColumns )
    AddColumns = list( set(AddColumns) - set( columnList ) )
    if len(AddColumns) != n:
        print('column already included!')
        return False
    return True

        
def PrintOvito( cordc, sfile, footer, attr_list=['x','y','z'] ):
    smat = cordc[attr_list]
    np.savetxt(sfile,[len(smat)],fmt='%s', footer='%s, %s'%(footer,str(attr_list)))
    np.savetxt(sfile,smat,fmt='%s')
    sfile.close()
    
def PltBitmap( value,
              xlabel = 'x', ylabel = 'y',
              xlim = (-0.5,0.5), ylim = (-0.5,0.5),
              frac = 1.0, #--- plot a patch
              zscore = True,
              title = 'cxy.png',
              colorbar=False,
              **kwargs
             ):
    fontsize = kwargs['fontsize'] if 'fontsize' in kwargs else 16
#    print('fontsize=',fontsize)
        
    val = value.copy()
    if 'scale' in kwargs and kwargs['scale'] == 'log':
        val = np.log10(val)
    #--- z-score
    if zscore:
#        print('mean=',np.mean(val))
#        print('std=',np.std(val))
        val -= np.mean(val)
        val /= np.std(val)
        val[val>2.0]=1.0
        val[val<-2.0]=-1.0
    vmin = kwargs['vmin'] if 'vmin' in kwargs else np.min(val)
    vmax = kwargs['vmax'] if 'vmax' in kwargs else np.max(val)
    #--- plot
    (mgrid,ngrid) = val.shape
    center = (ngrid/2,mgrid/2)
    #
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(labelsize=fontsize,which='both',axis='both', top=True, right=True)
    if 'labels' in kwargs:
        ax.axes.xaxis.set_visible(kwargs['labels'])
        ax.axes.yaxis.set_visible(kwargs['labels'])
    #
    pos = ax.imshow(val.real,cmap='bwr' if 'cmap' not in kwargs else kwargs['cmap'],
                     extent=(xlim[0],xlim[1],ylim[0],ylim[1]),origin='lower',
                    vmin=vmin, vmax=vmax,interpolation=kwargs['interpolation'] if 'interpolation' in kwargs else None )
    if colorbar:
        fig.colorbar( pos, fraction = 0.04)
    #
    ax.set_xlim(xlim[0]*frac,xlim[1]*frac)
    ax.set_ylim(ylim[0]*frac,ylim[1]*frac)
    #
        
    #
    PutMinorTicks(ax) #,MINUS=True)
    #
    if 'DrawFrame' in kwargs: 
        DrawFrame(ax, *kwargs['DrawFrame'])
    plt.savefig(title,dpi=2*75,bbox_inches='tight',pad_inches=0.0)
    plt.show()
    
    
def PltCrltn( value,
              xlabel = 'x', ylabel = 'y',
              xlim = (-0.5,0.5), ylim = (-0.5,0.5),
              frac = 1.0, #--- plot a patch
              zscore = True,
              fileName = 'cxy.png',
              dpi=75,
            ):
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    #
    val = value.copy()
    #--- zscore
    if zscore:
        val -= np.mean(val)
        val /= np.std(val)
        val[val>2.0]=1.0
        val[val<-2.0]=-1.0
    #
    (mgrid,ngrid) = val.shape
    center = (ngrid/2,mgrid/2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.axes.xaxis.set_visible(False) #--- remove labels
    ax.axes.yaxis.set_visible(False)

    pos = ax.imshow((CenterMatrix( val ).real),cmap='bwr',
                     extent=(xlim[0],xlim[1],ylim[0],ylim[1]),
                     #,vmin=-.01, vmax=.01
                    )
    ax.set_xlim(xlim[0]*frac,xlim[1]*frac)
    ax.set_ylim(ylim[0]*frac,ylim[1]*frac)

#    plt.colorbar( pos, fraction = 0.04)
    plt.savefig(fileName,dpi=dpi,bbox_inches='tight')
    plt.show()
    
def GetAutoCorrelation( val ):
    value  = val.copy()
    value -= np.mean( value )
    value /= np.std( value )

    ( nx, ny, nz ) =  value.shape
    n = nx * ny * nz
    vq = np.fft.fftn(value) #, axes=(0,1,2))
    vq_sq = np.abs(vq)**2

    v_real = np.fft.ifftn( vq_sq) / n #,  axes=(0,1,2) )
    return v_real

def CenterMatrix(a):
    ( mgrid, ngrid ) = a.shape
    return np.array([[a[i,j] for j in range(-int(ngrid/2),int(ngrid/2)+ngrid%2)] 
                              for i in range(-int(mgrid/2),int(mgrid/2)+mgrid%2)])

def Get_rc( xdata, ydata, cTOL ):
    try:
#        xc1 = xdata[np.abs(ydata)<cTOL][0] 
        xc1 = xdata[ydata<cTOL][0] 
    except:
        xc1 = np.nan
    try:    
        xc2 = xdata[ydata<0.0][0] #--- correlation length
    except:
        xc2 = np.nan
    try:
        xc = np.array([xc1,xc2])
        xc = xc[~np.isnan(xc)].min()
    except:
        xc = np.nan
    return xc

def PltCrltnFunc( crltn, 
                 xv,yv,
                 cTOL = 1.0e-2,
                 PLOT = True,
                 fileName='cxy.png',
                 title = '',
                 dpi = 60,
                 symmetry = True,
                ):
    (ny,nx,nz) = crltn.shape
    if PLOT:
        fig = plt.figure( figsize = (4,4))
        ax = fig.add_subplot(111)
        xstr = r'$r$'
        ystr = r'$c(r)$'
#         ax.set_xlabel(xstr,fontsize=16)
#         ax.set_ylabel(ystr,fontsize=16)
        ax.tick_params(labelsize=16)
        PutMinorTicks(ax)
        ax.tick_params(labelsize=20,which='both',axis='both', top=True, right=True)
    #
    val = crltn[:,:,0].copy() #--- only xy plane
    (m,n)=val.shape
    #--- along x 
    xv2 = xv[:,:,0].copy()
    dx = xv2[0,1] - xv2[0,0]
    #
    indx0 = 0 if symmetry else -int(n/2)
    rangee = np.arange(indx0,(int(n/2)+n%2))
    xdata0 = rangee * dx
    ydata0 = val[0,rangee]
    #
    xc = Get_rc( xdata0, ydata0, cTOL )
    #
    if PLOT:
        ax.plot( xdata0, ydata0,
                '-o',label=r'$x$',
                markersize=10,
                color='black',
                markerfacecolor='white',
                markeredgecolor=None,
                markevery=int(len(xdata0)/10),
               )       
    #--- along y 
    yv2 = yv[:,:,0].copy()
    dy = yv2[1,0] - yv2[0,0]
    #
    indx0 = 0 if symmetry else -int(m/2)
    rangee = np.arange(indx0,(int(m/2)+m%2))
    xdata = rangee * dy
    ydata = val[rangee,0]
    #
    yc = Get_rc( xdata, ydata, cTOL )
    #
    if PLOT:
        ax.plot( xdata, ydata,
                '-s', 
                color = 'red',
                label=r'$y$',
                markersize=10,
                markerfacecolor=None,
                markeredgecolor='black',
                markevery=int(len(xdata)/10),
               )
    #--- plot
    if PLOT:
        ax.legend(frameon=False, fontsize=20,handletextpad=.4,handlelength=1)
        ax.set_title( title )
        #ax.set_yscale('log')
        #
        ax.plot(ax.axis()[:2],[0.0,0.0],'-.',color='red')
#        ax.plot([dx,dx],ax.axis()[2:],'-.',color='black')
#        ax.plot((xc,xc),ax.axis()[2:],'-.r')
        #
        DrawFrame(ax, 0.2,0.09,0.15,0.06,0.04)
        #
        plt.savefig(fileName,dpi=2*75,bbox_inches='tight',pad_inches=0.0)
        plt.show()
    #
    return (xc, yc), (xdata0,ydata0), (xdata,ydata) 


def GetSlice2d( hist, err,
         xv, yv, zv,
         xtol = 2.5,
         z = 0.0):
############################################################
####### Get a 2D slice out off a 3D matrix
############################################################    
    dx=xtol #--- discretization length
    (ny,nx,nz) = hist.shape
    indices = np.all([np.abs(zv-z)>0.0,np.abs(zv-z)<dx],axis=0) #--- filtering based on the given range
    #--- binning in xy
    flist = hist[indices]
    rvect = np.c_[yv[indices],xv[indices]]
    rx,    bin_edges = np.histogramdd( rvect, bins = (ny, nx), weights = xv[indices] ) #--- \sum r_i
    ry,    bin_edges = np.histogramdd( rvect, bins = (ny, nx), weights = yv[indices]) #--- \sum r_i
    error,    bin_edges = np.histogramdd( rvect, bins = (ny, nx), weights = err[indices] ) #--- \sum r_i
    fmean, bin_edges = np.histogramdd( rvect, bins = (ny, nx), weights = flist ) #--- 3d histogram
    count, bin_edges = np.histogramdd( rvect, bins = (ny, nx) ) #--- n_i

    count[count==0] = 1
    rx /= count 
    ry /= count 
    fmean /= count
    error /= count
    return rx, ry, fmean, error

def GetSlice1d( hist, err,
         xv, yv,
         xtol = 2.5,
         **kwargs):
############################################################
####### Get a 1D slice out off a 2D matrix
############################################################    
    dx=xtol #--- discretization length
    (ny,nx) = hist.shape
    if 'y' in kwargs:
        y = kwargs['y']
        indices = np.all([np.abs(yv-y)>0.0,np.abs(yv-y)<dx],axis=0) #--- filtering based on the given range
        flist = hist[indices]
        rvect = xv[indices]
        rdist,    bin_edges = np.histogram( rvect, bins = nx, weights = xv[indices] ) #--- \sum r_i
        error,    bin_edges = np.histogram( rvect, bins = nx, weights = err[indices] ) #--- \sum r_i
        count, bin_edges = np.histogram( rvect, bins = nx ) #--- n_i
        fmean, bin_edges = np.histogram( rvect, bins = nx, weights = flist ) #--- 3d histogram
    elif 'x' in kwargs:
        x = kwargs['x']
        indices = np.all([np.abs(xv-x)>0.0,np.abs(xv-x)<dx],axis=0) #--- filtering based on the given range
        flist = hist[indices]
        rvect = yv[indices]
        rdist,    bin_edges = np.histogram( rvect, bins = ny, weights = yv[indices] ) #--- \sum r_i
        error,    bin_edges = np.histogram( rvect, bins = ny, weights = err[indices] ) #--- \sum r_i
        count, bin_edges = np.histogram( rvect, bins = ny ) #--- n_i
        fmean, bin_edges = np.histogram( rvect, bins = ny, weights = flist ) #--- 3d histogram
        
    #--- binning in xy

    count[count==0] = 1
    rdist /= count 
    fmean /= count
    error /= count
    return rdist, fmean, error


def PltCrltnFunc1d( crltn, err,
                 xv,
                 cTOL = 1.0e-2,
                 PLOT = True,
                 fileName='cxy.png',
                 title = '',
                 dpi = 60,
                 ylim=(-1.0,+1.0),
                ):
    if PLOT:
        fig = plt.figure( figsize = (4,4))
        ax = fig.add_subplot(111)
        xstr = r'$r$'
        ystr = r'$c(r)$'
        ax.set_xlabel(xstr,fontsize=16)
        ax.set_ylabel(ystr,fontsize=16)
        ax.tick_params(labelsize=16)
    #
    val = crltn.copy() #--- only xy plane
#    (m,n)=val.shape
    #--- along x 
    xv2 = xv.copy()
    dx = xv2[1] - xv2[0]
    #
    xdata = xv2 #np.arange(0,(n/2+n%2)) * dx
    ydata = val #[0,0:(n/2+n%2)]
    #
    xc = Get_rc( xdata, ydata, cTOL )
    #
    if PLOT:
        ax.errorbar( xdata, ydata,yerr = err, fmt='-o',label=r'$x$')       
    #--- plot
    if PLOT:
        ax.legend(frameon=False)
        ax.set_title( title )
        #ax.set_yscale('log')
        #
        ax.plot(ax.axis()[:2],[0.0,0.0],'-.',color='black')
#        ax.plot([dx,dx],ax.axis()[2:],'-.',color='black')
        ax.plot([xc,xc],ax.axis()[2:],'-.',color='black')
        #
#        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        #
        plt.savefig(fileName,dpi=dpi,bbox_inches='tight')
        plt.show()
    #
    return xc


def DrawFrame(ax, alpha_xl,alpha_xr,alpha_yb,alpha_yt,linewidth,LOG_X=None,LOG_Y=None):
    [xlo,xhi,ylo,yhi] = ax.axis()
    if LOG_X:
        [xlo,xhi,junk,junk] = np.log10(ax.axis())
    if LOG_Y:
        [junk,junk,ylo,yhi] = np.log10(ax.axis())
    lx = xhi - xlo
    ly = yhi - ylo
    xy = [xlo - alpha_xl * lx, ylo - alpha_yb * ly]
    height = ly*(1+alpha_yb+alpha_yt)
    width = lx*(1+alpha_xl+alpha_xr)
    xy_end=[xy[0]+width,xy[1]+height]
    if LOG_X:
        xy[0] = 10 ** xy[0]
        xy_end[0] = 10 ** xy_end[0]
    if LOG_Y:
        xy[1] = 10 ** xy[1]
        xy_end[1] = 10 ** xy_end[1]
    ax.add_patch( patches.Rectangle(xy=xy, width=xy_end[0]-xy[0], 
                                    height=xy_end[1]-xy[1], linewidth=linewidth,
                                    clip_on=False,facecolor=None,edgecolor='black',fill=None) ) 
    
def MultipleFrames( path='', title='', irun = 0, nmax = 10000 ):
    i=0
    append = False
    while i < nmax:
        try:
            sarr0 = np.c_[np.loadtxt('%s%i/Run%s/%s'%(path,i,irun,title))].T
            #        print i,sarr0
            if not append:
                sarr = sarr0.copy()
                append = True
            else:
                sarr = np.concatenate((sarr,sarr0),axis=0)
        except:
            i+=1
    #        traceback.print_exc()
            continue
        i+=1
    return sarr

def MultipleFrames2nd( path='', title='', nrun = 0, ncols=3 ):
    i=0
    append = False
#    pdb.set_trace()
    while i < nrun:
        sarr0 = (np.ones(ncols)*np.nan).reshape((1, ncols))
#        print(sarr0.shape)
        try:
            sarr0 = np.c_[np.loadtxt('%s/Run%s/%s'%(path,i,title))].T
#            print(sarr0.shape)
        except:
#            traceback.print_exc()
            pass
        if not append:
            sarr = sarr0.copy()
            append = True
#            print(i,sarr0)
        else:
#            print(i,sarr0)
            sarr = np.concatenate((sarr,sarr0),axis=0)
        i+=1
    return sarr

def AvalancheSize(strain, Virial):
    d={'t_end':[],'duration':[],'ds':[]}
    kount = 0
    duration = 0
    ds = 0.0
#    counter = 0
    try:
        dtt = strain[1]-strain[0] #TimeSeries[0]['Time'].iloc[1]-TimeSeries[0]['Time'].iloc[0]
        for items,sbulk in list(zip(strain, Virial)): #TimeSeries[isamp].itertuples():
#            sbulk = items[2]
            t = items #items[1]
            #--- ens. average
            if kount == 0:
                a = sbulk
                ta = t
            elif kount == 1:
                b = sbulk
                tb = t
            else:
                c = sbulk
                tc = t
            if kount < 2:
                kount += 1
                init0 = kount + 1
                continue
            sdr = 0.5 * ( c - a ); #--- derivative
#                print(t, sdr)
#            if sdr < 0.0 and counter == 0:
#                continue
            if sdr > 0.0: #--- elastic loading
                init0 = kount + 1; #--- initialize init
            else: #--- avalanche starts!
#                    print(t, sdr)
                ds += sdr #--- increment avalanche size by the stress derivative
                duration += 1 #--- increment duration
            if init0 - kount == 1 and duration != 0: #--- avalanche ends!
                print(duration,ds) #tc-duration*(tb-ta),tc,duration
    #			ax.plot([tc-duration*(tb-ta),tc], [0.0,0.0]
    #                    ,'o')
                d['t_end'].append(tc) 
                d['duration'].append(duration*dtt) 
                d['ds'].append(-ds)
                ds = 0.0 #--- initialize 
                duration = 0
#                counter+=1
    #		if counter == 1:
    #			break
            a = b #--- new assignments
            b = c;
            ta = tb #--- new assignments
            tb = tc;
            kount += 1
    except:
#        traceback.print_exc()
        pass
    if duration != 0: #--- one single avalanche!
        d['t_end'].append(tc) 
        d['duration'].append(duration*dtt) 
        d['ds'].append(-ds)
#    print(duration,ds)
#fig.show()
    df=pd.DataFrame(d)
    df=df[df['ds']!=0.0]
    
    return df

def GetPDF(slist, n_per_decade=4, ACCUM = None, linscale = None, density=True):
    if not linscale:
        xlo = np.floor(np.log10(np.min(slist)))
        xhi = np.ceil(np.log10(np.max(slist)))
        bins = np.logspace(xlo,xhi,int(xhi-xlo)*n_per_decade)    
    else:
        xlo = np.min(slist)
        xhi = np.max(slist)
        bins = np.linspace(xlo,xhi,n_per_decade)
        
    hist, edges = np.histogram(slist,bins=bins,density=density)
    count, edges = np.histogram(slist,bins=bins)
    
    
    if ACCUM:
        return np.cumsum((edges[1:]-edges[:-1])*hist), edges
    
    
    hist = hist[count>1]
    edges = edges[:-1][count>1]
    count = count[count>1]
    
    return  hist, edges, hist / count**0.5


# create a definition for the short hyphen
#matplotlib.rcParams["text.latex.preamble"]+= r'\mathchardef\mhyphen="2D'#.append(r'\mathchardef\mhyphen="2D')

class MyLogFormatter(matplotlib.ticker.LogFormatterMathtext):
    def __call__(self, x, pos=None):
#        pass
        # call the original LogFormatter
        rv = matplotlib.ticker.LogFormatterMathtext.__call__(self, x, pos)

        # check if we really use TeX
        if matplotlib.rcParams["text.usetex"]:
            # if we have the string ^{- there is a negative exponent
            # where the minus sign is replaced by the short hyphen
            rv = re.sub(r'\^\{-', r'^{\mhyphen', rv)

        return rv
    
def PutMinorTicks(ax, LOGY=None,LOGX=None, MINUS=False):
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    if LOGY:
        #--- add major yticks
        ymin=np.ceil(np.log10(ax.axis()[2]))
        ymax=np.floor(np.log10(ax.axis()[3]))
        nbin = ymax - ymin
        ax.set_yticks(np.logspace(ymin,ymax,int(nbin)+1))
        #--- put minor bins y
        locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),numticks=12)
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if LOGX:
        #--- add major yticks
        ymin=np.ceil(np.log10(ax.axis()[0]))
        ymax=np.floor(np.log10(ax.axis()[1]))
        nbin = ymax - ymin
        ax.set_xticks(np.logspace(ymin,ymax,int(nbin)+1))
        #--- put minor bins y
        locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),numticks=12)
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    if MINUS:
        ax.xaxis.set_major_formatter(MyLogFormatter()) #--- minus sign too long
        ax.yaxis.set_major_formatter(MyLogFormatter()) #--- minus sign too long

        
def makeTickLabels(ax, xdata, ydata, **kargs):
    ylo = kargs['ylo'] if 'ylo' in kargs else 10**np.floor(np.log10(np.min(ydata)))
    yhi = kargs['yhi'] if 'yhi' in kargs else 10**np.ceil(np.log10(np.max(ydata)))
    xlo = kargs['xlo'] if 'xlo' in kargs else 10**np.floor(np.log10(np.min(xdata)))
    xhi = kargs['xhi'] if 'xhi' in kargs else 10**np.ceil(np.log10(np.max(xdata)))
    center = kargs['center'] if 'center' in kargs else True
    MINUS = kargs['MINUS'] if 'MINUS' in kargs else True
    xc = 0.5*np.log10(xhi*xlo)
    yc = 0.5*np.log10(yhi*ylo)
    
    dx = np.log10(xhi/xlo)*0.5
    dy = np.log10(yhi/ylo)*0.5

    if center:
        dx = dy = np.max([dx,dy])
    
    ax.axis(10**np.array([xc-dx,xc+dx,yc-dy,yc+dy]))
    ax.loglog()
    
    ##--- add major xticks
    xmin=np.ceil(np.log10(ax.axis()[0]))
    xmax=np.floor(np.log10(ax.axis()[1]))
    nbin = xmax - xmin
    ax.set_xticks(np.logspace(xmin,xmax,int(nbin)+1))
    
    #--- add major yticks
    ymin=np.ceil(np.log10(ax.axis()[2]))
    ymax=np.floor(np.log10(ax.axis()[3]))
    nbin = ymax - ymin
    ax.set_yticks(np.logspace(ymin,ymax,int(nbin)+1))
    
    #--- put minor bins
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),numticks=12)
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if MINUS:
        ax.xaxis.set_major_formatter(MyLogFormatter()) #--- minus sign too long

    #--- put minor bins y
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),numticks=12)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    if MINUS:
        ax.yaxis.set_major_formatter(MyLogFormatter()) #--- minus sign too long

        
    ax.tick_params(axis='y',left=True, right=True,which='both')
    ax.tick_params(axis='x',bottom=True, top=True,which='both')
    
def GetBinnedAverage( a, y, **kwargs ):
    n=len(a)
    if 'nbins_per_decade' in kwargs:
        nbins = kwargs['nbins_per_decade'] * int( ( np.ceil(np.max(a))-np.floor(np.min(a)) ) )
    if 'nbins' in kwargs:
        nbins = kwargs['nbins']
    
    ysum = np.histogram(a, bins=nbins, weights=y)[0]
    ysum_sq = np.histogram(a, bins=nbins, weights=y*y)[0]
    xsum = np.histogram(a, bins=nbins, weights=a)[0]
    xsum_sq = np.histogram(a, bins=nbins, weights=a*a)[0]
    xcount = np.histogram(a, bins=nbins)[0]

#    print xsum
#    print xcount
#    assert not np.any(xcount==0)
    #--- remove zero
    xsum = xsum[xcount!=0]
    xsum_sq = xsum_sq[xcount!=0]
    ysum = ysum[xcount!=0] 
    ysum_sq = ysum_sq[xcount!=0] 
    xcount = xcount[xcount!=0]

    xmean = xsum/xcount
    ymean = ysum/xcount
    
    ystd2=ysum_sq/xcount-ymean*ymean
    xstd2=xsum_sq/xcount-xmean*xmean
    
    return xmean, ymean, 2*(xstd2 / xcount) ** 0.5, 2*(ystd2 / xcount) ** 0.5

def Intrp( d2min, box0, attr, Plot = None, title = 'test.png',**kwargs ):
    #--- mean dist between atoms 
    natoms = len( d2min.x ) 
    CellVectorOrtho, VectorNorm = lp.GetOrthogonalBasis( box0.CellVector )
    volume = np.linalg.det( CellVectorOrtho )
    dmean = 0.5*( volume / natoms ) ** (1.0/3.0) 


    #--- grid tiling mapped box with original size
    #--- values are interpolated onto this grid
    (xlin, ylin, zlin), (xv, yv, zv) = lp.GetCubicGrid( box0.CellOrigin, 
                                                     box0.CellVector, 
                                                     dmean,
                                                     margin = 0.0 * dmean )
    xi = np.array(list(zip(xv.flatten(), yv.flatten(), zv.flatten())))

    #--- expand the original box
        #--- map to square box
    mapp = lp.Map( d2min, box0 ) 
    mapp.ChangeBasis()
    mapp.Set( d2min ) #--- atoms: add mapped xyz

    cptmp = lp.Copy(d2min, box0) #--- important: must be reference frame!!
    cptmp.Expand( epsilon = 0.2, mode = 'isotropic' )
    d2exp = cptmp.Get()

    points = np.c_[d2exp.xm,d2exp.ym,d2exp.zm] #--- unstructured points
    values = np.array(d2exp[attr]) #(np.array(d2exp.C66)+np.array(d2exp.C55)+np.array(d2exp.C44))/3.0 #np.c_[-(np.array(d2exp.sxx)+np.array(d2exp.syy)+np.array(d2exp.szz))/3.0/np.array(d2exp.AtomicVolume)] #--- corresponding values

    grid_z = scp_int.griddata(points, values, xi, method='linear')
    assert not np.any(np.isnan(grid_z.flatten())), 'increase ev!'

    #--- make an object
    d2intrp = lp.Atoms(**pd.DataFrame(np.c_[xi,grid_z],columns=['x','y','z',attr]).to_dict(orient='series'))

    #--- reshape value
    nx,ny,nz = len(xlin), len(ylin),len(zlin) 
    value = np.c_[d2intrp[attr]].reshape(((ny,nx,nz)))

    CellVectorOrtho, VectorNorm = lp.GetOrthogonalBasis( box0.CellVector ) #--- box length

    #--- xy plane
    #--- 2d slice
    nzl=[0]
    val = Get2dSlice( value, zlin, 
                        zlin[-1], nzll=nzl  )
    nzz=nzl[0]
        #
    #--- square bitmap
    lx=np.min([xlin[-1]-xlin[0],ylin[-1]-ylin[0]])
    xc = 0.5*(xlin[-1]+xlin[0])
    yc = 0.5*(ylin[-1]+ylin[0])

    if Plot:
        PltBitmap(val, 
#                  xlim=VectorNorm[0]*np.array([-0.5,0.5]),ylim=VectorNorm[1]*np.array([-0.5,0.5]),
                  xlim=np.array([xc-0.5*lx,xc+0.5*lx]),ylim=np.array([yc-0.5*lx,yc+0.5*lx]),
                  frac = 1.0, #--- plot a patch
                  title = title,
                  **kwargs
                )

#    return (xlin, ylin, zlin), (xv[:,:,nzz], yv[:,:,nzz], zv[:,:,nzz]), d2intrp
    return (xlin, ylin, zlin), (xv, yv, zv), d2intrp


def IntrpScatter( d2min, box0, attr, vor, **kwargs ):
    #--- mean dist between atoms 
    natoms = len( d2min.x ) 
    CellVectorOrtho, VectorNorm = lp.GetOrthogonalBasis( box0.CellVector )
    volume = np.linalg.det( CellVectorOrtho )
    dmean = 0.5*( volume / natoms ) ** (1.0/3.0) 


    #--- grid tiling mapped box with original size
    #--- values are interpolated onto this grid
    (xlin, ylin, zlin), (xv, yv, zv) = lp.GetCubicGrid( box0.CellOrigin, 
                                                     box0.CellVector, 
                                                     dmean,
                                                     margin = 0.0 * dmean )
    xi = np.array(list(zip(xv.flatten(), yv.flatten(), zv.flatten())))

    #--- expand the original box
        #--- map to square box
    mapp = lp.Map( d2min, box0 ) 
    mapp.ChangeBasis()
    mapp.Set( d2min ) #--- atoms: add mapped xyz

    cptmp = lp.Copy(d2min, box0) #--- important: must be reference frame!!
    cptmp.Expand( epsilon = 0.2, mode = 'isotropic' )
    d2exp = cptmp.Get()

    points = np.c_[d2exp.xm,d2exp.ym,d2exp.zm] #--- unstructured points
    values = np.array(d2exp[attr]) #(np.array(d2exp.C66)+np.array(d2exp.C55)+np.array(d2exp.C44))/3.0 #np.c_[-(np.array(d2exp.sxx)+np.array(d2exp.syy)+np.array(d2exp.szz))/3.0/np.array(d2exp.AtomicVolume)] #--- corresponding values

    grid_z = scp_int.griddata(points, values, xi, method='linear')
    assert not np.any(np.isnan(grid_z.flatten())), 'increase ev!'

    #--- make an object
    d2intrp = lp.Atoms(**pd.DataFrame(np.c_[xi,grid_z],columns=['x','y','z',attr]).to_dict(orient='series'))

    if 'Plot' in kwargs and kwargs['Plot']:


        
    #---    
    #--- background map
    #---
    #--- reshape value
        nx,ny,nz = len(xlin), len(ylin),len(zlin) 
        value = np.c_[d2intrp[attr]].reshape(((ny,nx,nz)))
        #
        CellVectorOrtho, VectorNorm = lp.GetOrthogonalBasis( box0.CellVector ) #--- box length
        
        
        #--- get xy plane
        #--- 2d slice
        val = Get2dSlice( value, zlin, 
                            zlin[-1],
                            #0.5*(zlin[0]+zlin[-1]),
                        )
        
        #---
        #--- get scatterd data
        #---
        mapp = lp.Map( vor, box0 ) 
        mapp.ChangeBasis()
        mapp.Set( vor ) #--- atoms: add mapped xyz
        #--- filter
        indices = GetFullIcosahedra( vor )

            #--- plot scattered data
        xs = np.array(vor.xm)[indices]
        ys = np.array(vor.ym)[indices]
        zs = np.array(vor.zm)[indices]
        xyzv = Get2dSliceScattered( xs, ys, zs, 
                                    zlin, 
                                   zlin[-1],
                                   #0.5*(zlin[0]+zlin[-1])
                                  )
#        print('xmin=',xyzv[:,0].min())
#        print('xlin=',np.min(xlin))
    
        #
#        dz *= 2.0
#        zlo=(zz-dz-zlin[0]) % lz + zlin[0]
#        zhi=(zz+dz-zlin[0]) % lz + zlin[0]
#        print(zlo,zc,zhi)
#        pdb.set_trace()
#        indxx= np.all([vor.zm > zlo, 
#                       vor.zm < zhi
#                      ],
#                        axis=0)
#        if zhi < zlo: #--- periodic boundaries
#            indxx= np.any([vor.zm > zlo, 
#               vor.zm < zhi
#              ],
#                axis=0)
#        indxx = np.all([indxx,indices],axis=0)
#        xyzv = np.c_[np.array(vor.xm)[indxx],np.array(vor.ym)[indxx]]

        #--- square bitmap
        lx=np.min([xlin[-1]-xlin[0],ylin[-1]-ylin[0]])
        xc = 0.5*(xlin[-1]+xlin[0])
        yc = 0.5*(ylin[-1]+ylin[0])

        #--- bitmap
        PltBitmapWithScatter(val, 
              xyzv,
#                  xlabel = r'$x$(\r{A})', ylabel = r'$y$(\r{A})',
                  xlim=np.array([xc-0.5*lx,xc+0.5*lx]),ylim=np.array([yc-0.5*lx,yc+0.5*lx]),
#                  xlim=np.array([xlin.min(),xlin.max()]),ylim=np.array([ylin.min(),ylin.max()]),
              **kwargs
#                     DrawFrame=[0.2,0.09,0.15,0.06,0.04],
              )
            
            
#         PltBitmap(val, 
#                   xlabel = r'$x$(\r{A})', ylabel = r'$y$(\r{A})',
#                   xlim=VectorNorm[0]*np.array([-0.5,0.5]),ylim=VectorNorm[1]*np.array([-0.5,0.5]),
#                   zscore = True,
#                   frac = 1.0, #--- plot a patch
#                   title = 'd2min.png',
#                   colorbar=True,
#                 )

    return (xlin, ylin, zlin), (xv, yv, zv), d2intrp


def GetAverage( X, y, nbins ):
    assert X.shape == y.shape
    
    hist, edges = np.histogram(X.flatten(), bins=nbins, weights = y.flatten())
    sum_sq, edges = np.histogram(X.flatten(), bins=nbins, weights = y.flatten()**2)
    count, edges = np.histogram(X.flatten(), bins=nbins)

    hist = hist[count>0]
    sum_sq = sum_sq[count>0]
    edges = edges[:-1][count>0]
    count = count[count>0]

    hist /= count
    sum_sq /= count
    xstd = (sum_sq-hist**2)**0.5

    return edges,hist, xstd


def sroDensityDiscrete( vor, box,
                AtomicRadius,     
                 **kwargs,
               
                ):

        #--- filter
        indices = GetFullIcosahedra( vor )        
        #--- output dump file (apply filtering), run ovitos
        atomFiltrd = lp.Atoms(**pd.DataFrame(vor.__dict__)[indices].to_dict(orient='series')) #--- filter atoms
        #
        wdf = lp.WriteDumpFile(atomFiltrd,box) #--- need high-precision data
        wdf.Write('junk.xyz')
        #--- load to ovitos
        rad1=0.0#AtomicRadius[1]
        rad2=0.0#AtomicRadius[2]
        rad3=0.0#AtomicRadius[3]
        rad4=0.0#AtomicRadius[4]
        rad5=0.0#AtomicRadius[5]
        os.system("ovitos OvitosCna.py \'junk.xyz\' \'VoronoiFiltrd.xyz\' 1 3 %s %s %s %s %s"%(rad1,rad2,rad3,rad4,rad5)) 
    
        #--- read from ovito output
        lmpDataFiltrd = lp.ReadDumpFile( 'VoronoiFiltrd.xyz' )
        lmpDataFiltrd.GetCords( ncount = sys.maxsize )
        #
        vorFiltrd = lp.Atoms( **lmpDataFiltrd.coord_atoms_broken[0].to_dict(orient='series') )
        #--- filter d2min
 #       d2minFiltrd[itimee] = lp.Atoms(**pd.DataFrame(d2min.__dict__)[indices].to_dict(orient='series')) #--- filter atoms

        #--- density
        x = 1.0/np.array(vorFiltrd.AtomicVolume)

        #--- reformat the array
        natoms = len(vor.x)
        rho = np.zeros(natoms)
        rho[indices] = x        
        return rho
        
