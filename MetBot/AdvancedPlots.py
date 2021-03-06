# AdvancedPlots.py
import SynopticAnatomy as sy
import numpy as np
import matplotlib.pyplot as plt
import EventStats as stats
#import SynopticClassify as sc
import SynopticPlot as syp
import RainStats as rs
import mytools as my
import MetBlobs as blb
from time import sleep as sleep

### SUBPLOT FUNCTIONS
def spatiofreq2_season(s,lat,lon,yrs,eventkeys,figno=1,season='coreseason',\
     key='noaa-olr-0-all',flagonly=False,file_suffix='test',savefig=False,\
     fontdict=False):
    '''spatiofreq2_season(s,lat,lon,yrs,figno=1,season='coreseason',\
                          flagonly=False,file_suffix='test')
    Produces subplots of cloud-band gridpoint count by month
    '''
    if not fontdict:fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbkl=key.split('-')
    if mbkl[0]=='noaa':dclim=(1,9,1)
    if mbkl[0]=='hadam3p':dclim=(1,11,1)

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season
    m, f = blb.SAfrBasemap(lat[4:-7],lon[3:-2],drawstuff=True,prj='cyl',\
           fno=figno,rsltn='l',fontdict=fd)
    if len(mns)==12:plt.close();plt.figure(figsize=[12,10])
    elif len(mns)==6:plt.close();plt.figure(figsize=[13,8])
    cnt=1
    msklist=[]
    for mn in mns:
        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)
        if flagonly:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,\
                    clim=(1.,9.,1.),month=mn,flagonly=True,fontdict=fd)
        else:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,clim=dclim,\
                                      month=mn,flagonly=False,fontdict=fd)
        if cnt%2==0:
            syp.redrawmap(m,lns=True,resol='verylow',parallel=False,fontdict=fd)
        else:
            syp.redrawmap(m,lns=True,resol='verylow',fontdict=fd)
        cnt+=1
        msklist.append(allmask)
        my.xtickfonts();my.ytickfonts()
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(stats.figdir+file_suffix+'_flagonly.png',dpi=150)
        else: plt.savefig(stats.figdir+file_suffix+'.png',dpi=150)

    return msklist

def spatiofreq2_seasonanoms(s,lat,lon,yrs,eventkeys,msklist,figno=1,\
    season='coreseason',key='noaa-olr-0-all',flagonly=False,\
    file_suffix='test',savefig=False,fontdict=False):
    '''spatiofreq2_season(s,lat,lon,yrs,figno=1,season='coreseason',\
                          flagonly=False,file_suffix='test')
    Produces subplots of cloud-band gridpoint count by month
    '''
    if not fontdict: fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbkl=key.split('-')
    if mbkl[0]=='noaa':dclim=(-3.,3,.5)
    if mbkl[0]=='hadam3p':dclim=(-5.,5,1.)

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season
    m, f = blb.SAfrBasemap(lat[4:-7],lon[3:-2],drawstuff=True,prj='cyl',\
                           fno=figno,rsltn='l',fontdict=fd)
    if len(mns)==12:plt.close();plt.figure(figsize=[12,10])
    elif len(mns)==6:plt.close();plt.figure(figsize=[13,8])
    cnt=1
    for mn in mns:
        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)
        if flagonly:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,clim=dclim,\
            meanmask=msklist[cnt-1],month=mn,flagonly=True,fontdict=fd)
        else:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,fontdict=fd,\
            meanmask=msklist[cnt-1],clim=dclim,month=mn,flagonly=False)
        if cnt%2==0:
            syp.redrawmap(m,lns=True,resol='verylow',parallel=False,\
                          fontdict=fd,clr='0.5')
        else:
            syp.redrawmap(m,lns=True,resol='verylow',fontdict=fd,clr='0.5')
        cnt+=1
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(stats.figdir+file_suffix+'_cbcount_flagonly.png',\
                        dpi=150)
        else:
            plt.savefig(stats.figdir+file_suffix+'_cbcount.png',dpi=150)
