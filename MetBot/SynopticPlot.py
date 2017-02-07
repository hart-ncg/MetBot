'''SynopticPlot.py: A module of the MetBot package

This module contains a collection of functions to automate plotting of
1: events described by SynopticAnatomy.Events objects
2: composite events computed from some of these events'''

import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import mytools as my
import mynetcdf as mync
import MetBlobs as blb
import SynopticAnatomy as sy
import EventStats as stats
import matplotlib.path as mpath
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from time import sleep as sleep
import threading
typlt=type(plt)
tymp=bm.Basemap
# PLOT PROPERTIES
lineprops = {0:{'trk': "'k-',lw=2.0" , 'blob': "color='0.1',ls='--',lw='2.0'"},
             1:{'trk': "'r-',lw=2.0" , 'blob': "color='0.5',ls='--',lw='2.0'"},
             2:{'trk': "'c-',lw=2.0" , 'blob': "color='0.8',ls='--',lw='2.0'"},
             3:{'trk': "'m-',lw=2.0" , 'blob': "color='k',ls='--',lw='3.0'"},
             4:{'trk': "'k-',lw=2.0" , 'blob': "color='k',ls='-.',lw='3.0'"},}

lineprops = {0:{'trk': "'k-',lw=2.0" , 'blob': "color='0.4',ls='--',lw='3.0'"},
             1:{'trk': "'r-',lw=2.0" , 'blob': "color='m',ls='--',lw='3.0'"},
             2:{'trk': "'c-',lw=2.0" , 'blob': "color='0.8',ls='--',lw='2.0'"},
             3:{'trk': "'m-',lw=2.0" , 'blob': "color='k',ls='--',lw='3.0'"},
             4:{'trk': "'k-',lw=2.0" , 'blob': "color='k',ls='-.',lw='3.0'"},}
gcols = ['0.1', '0.5', '0.8'] # Grey scale intensties
cols =  ['r','b','k']
tcols =  ['c','m','0.6']
# DATA LOADING
def openolr(reftime=False,sub='SA'):
    '''Directory tree and file specific function for opening gridded olr data used in MetBot modules'''
    hrwindow=49
    v="noaa-olr-0-0";dset, varstr, lev, drv = v.split('-')
    dsrc="/home/neil/data/olr/"
    dsrc="/home/neil/sogehome/data/olr/"
    olr, time, lat, lon, dtime = my.mync.openolr(dsrc+'olr.day.mean.nc','olr',subs='SA')
    time=time-15769752.0;olr=olr[1675:,:,:];time=time[1675:];dtime=dtime[1675:]
    if isinstance(reftime, np.ndarray): exec("ixt, [time, "+varstr+", dtime] = my.ixtwindow(reftime,time,hrwindow,time,"+varstr+",dtime)")
    return (olr,time,lat,lon,dtime)

def openmetvars(varlist,reftime=False,sub='SA',daily=False,\
                datadir="/home/neil/data/ncep2/6hourly/"):
    '''Directory tree and file specific function for opening gridded reanalysis data used in MetBot modules'''
    hrwindow=49
    varstrlist=[]
    for v in varlist:
        sbst=sub;dset, varstr, levsel, deriv = v.split('-')
        dsrc=datadir+varstr+"/"
        dsrc=datadir
        ncfile = dsrc+"%s.%s_%s.nc" %(varstr, sbst, levsel)
        if daily: ncfile = dsrc+"%s.%s_%s.daily.nc" %(varstr, sbst, levsel)
        exec(varstr+levsel+", time, lat, lon, dtime = my.mync.open"+dset+"(ncfile,\'"+varstr+"\')")
        if isinstance(reftime, np.ndarray): exec("ixt, [time, "+varstr+levsel+", dtime] = my.ixtwindow(reftime,time,hrwindow,time,"+varstr+levsel+",dtime)")
        varstrlist.append(varstr+levsel)
    exec('returnall=('+ ",".join(varstrlist)+',time,lat,lon,dtime)')
    return returnall

def openmetvars_flavours(varlist,enskeys,reftime=False,sub='SA'):
    '''Directory tree and file specific function for opening gridded reanalysis data used in MetBot modules'''
    hrwindow=49
    varstrlist=[]
    for v in varlist:
        sbst=sub;dset, varstr, levsel, deriv, exp = v.split('-')
        dset="ncep2"
        dsrc="/home/neil/data/flavours/"+exp+"/"
        #ncfileout = dsrc+"%s.%s.%s_%s.nc" %(exp,varstr, sbst, levsel)
        nclist=[dsrc+"%s.%s.%s_%s.nc" %(exp+k,varstr, sbst, levsel) for k in enskeys]
        exec(varstr+levsel+", time, lat, lon, lev, dtime = my.OpenMultipleNC(nclist,varstr,dset)")
        if isinstance(reftime, np.ndarray): exec("ixt, [time, "+varstr+levsel+", dtime] = my.ixtwindow(reftime,time,hrwindow,time,"+varstr+levsel+",dtime)")
        varstrlist.append(varstr+levsel)
    exec('returnall=('+ ",".join(varstrlist)+',time,lat,lon,dtime)')
    return returnall

# DATA PREPARATION
def mapxy(m,x,y,closeloop=False):
    '''use gcpoints to return values to plot things on the map
    This is only really needed for using '''
    lats,lons=np.array(()),np.array(())
    for i in xrange(len(x)-1):
        lt,ln = m.gcpoints(x[i],y[i],x[i+1],y[i+1],10)
        lats=np.append(lats,lt)
        lons=np.append(lons,ln)
    if closeloop:
        lt,ln = m.gcpoints(x[i],y[i],x[0],y[0],10)
    return lats, lons

def maskeddata(e,i,data,datatime,datakey):
    '''mdata = maskeddata(e,i,data,datatime,datakey)
    Returns a lat x lon array ready for plotting
    Usage: e is the event currently working with
           i is index of the time slice in the array
           data is a time x lat x lon array
           datatime is the time coordinates of data
           datakey is the mbs key of the variable
           tstep is temporal resolution of data
                 used because 24hr time res can be problem when using
                 24hr OLR with 6hr ncep2 data'''
    if isinstance(e.trkarrs[datakey],dict):
        print 'No track in %s for event %d at position %d' %(datakey,e.trkkey,i)
        data=np.zeros(data[0,:,:].shape,dtype=np.float32);data[:]=np.nan
        msk = np.ones(data.shape,dtype=np.bool)
        mdata = np.ma.MaskedArray(data,mask=msk)
        return mdata
    t = e.trkarrstime[datakey][i]
    #if tres==24: t = t - (t % 24)
    #if i<3 and tres==24 and (t % 24) > 0: t+=24
    #print datakey, t%24, my.ixtwindow([t],datatime,0)
    ix = my.ixtwindow([t],datatime,0)[0]
    mdata = np.ma.MaskedArray(data[ix,:,:],mask=e.masks[datakey][i,:,:])
    return mdata

def rainmaskeddata(e,i,data,datatime,datakey):
    '''rmdata = rainmaskeddata(e,i,data,datatime,datakey)
    Returns a masked station x 1 array ready for plotting
    Usage: e is the event currently working with
           i is index of the time slice in the array
           data is a station x time
           datatime is the time coordinates of data
           datakey is the mbs key of the variable'''
    t = e.trkarrstime[datakey][i]
    if t > datatime[-1]:
        print 'No more raindata to compute'
        return False
    ix = my.ixtwindow([t],datatime,0)[0]
    rmdata = np.ma.MaskedArray(data[:,ix],mask=e.rainmasks[datakey][i,:])
    return rmdata

def daylabel(e,etime):
    '''Welcome to the ugliest piece of code I've ever written...don't ask questions
    it just works OK! Siesa

    It returns a label for each time step in an event
    which describes the relative position of that step to
    flagged days'''
    masktimes = e.trkarrstime[e.refkey]
    # to ensure trkarrstime 6hr has valid values
    if isinstance(e.trkarrs[e.mbskeys[3]],np.ndarray):time6hr = e.trkarrstime[e.mbskeys[3]]
    elif isinstance(e.trkarrs[e.mbskeys[4]],np.ndarray):time6hr = e.trkarrstime[e.mbskeys[4]]
    elif isinstance(e.trkarrs[e.mbskeys[5]],np.ndarray):time6hr = e.trkarrstime[e.mbskeys[5]]
    new = np.zeros(masktimes.shape)
    dlabel = np.zeros(masktimes.shape)
    ddiff = np.zeros(masktimes.shape)
    dlabel[:]=-1
    for tf in e.flagtimes:
        new = np.where((masktimes==tf) & (new==0),tf,new)
    tdiff2=np.append(new[1:]-new[:-1],0)
    istrts = np.where(tdiff2>0)[0]
    dcnt=0.0
    for ist in istrts:
        dlabel[ist+1:ist+5] = dcnt
        dcnt+=1
    iscnt=0
    for ixt in xrange(len(time6hr)):
        t6 = time6hr[ixt]
        if dlabel[ixt] != -1: iscnt =  dlabel[ixt]
        t6diff = t6 - masktimes[istrts[iscnt]+1]
        ddiff[ixt] = t6diff
    ddiff = np.where(ddiff>18,ddiff-18,ddiff)
    flarr = np.int32(np.vstack((dlabel,ddiff)).T)
    lablist=[]
    for dlab, dd in flarr:
        if dlab!=-1: labl = 'D'+str(dlab+1)+'_'+str(dd)
        else:
            if dd<0: labl = str(dd)
            else: labl = '+'+str(dd)
        lablist.append(labl)

    if isinstance(etime,tuple):
        emin, emax = etime
        ixlt = np.where((ddiff>=emin) & (ddiff<=emax))[0]
        #ixlab = np.where(dlabel!=-1)
        #ixlist = np.unique(np.append(ixlt,ixlab))
        ixlist = ixlt
        ixlist = list(ixlist)
    elif isinstance(etime,list):
        ixlist=[]
        for flagd in etime:
            ixtrue = np.where( dlabel == flagd)[0]
            ixlist.extend(list(ixtrue))

    return np.asarray(lablist), np.asarray(ixlist), flarr[ixlist,:]

# PLOTTING FUNCTIONS
def redrawmap(m,lns=False,resol='low',parallel=True,meridian=True,\
    fontdict=False,clr='k'):
    '''Redraw the basemap features after calling plt.clf()'''
    if not fontdict: fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    m.drawcountries(color=clr)
    m.drawcoastlines(color=clr)
    if lns and resol=='low':
        delon = 10.; meridians = np.arange(10.,360.,delon)
        delat = 5.; circles = np.arange(0.,90.+delat,delat).tolist()+\
                              np.arange(-delat,-90.-delat,-delat).tolist()
    elif lns and resol=='hi':
        delon = 2.; meridians = np.arange(10.,360.,delon)
        delat = 2.; circles = np.arange(0.,90.+delat,delat).tolist()+\
                              np.arange(-delat,-90.-delat,-delat).tolist()
    elif lns and resol=='verylow':
        delon = 20.; meridians = np.arange(10.,360.,delon)
        delat = 10.; circles = np.arange(0.,90.+delat,delat).tolist()+\
                               np.arange(-delat,-90.-delat,-delat).tolist()
    if parallel:
        m.drawparallels(circles,linewidth='0.1',labels=[1,0,0,0],fontdict=fd)
    if meridian:
        m.drawmeridians(meridians,linewidth='0.1',labels=[0,0,0,1],fontdict=fd)

def addpatch(verts):
    '''Helper to add patch defined by poly with xverts, yverts'''
    #verts=np.hstack((xverts[:,np.newaxis],yverts[:,np.newaxis]))
    Path = mpath.Path
    moves = ','.join(['Path.LINETO' for i in xrange(verts.shape[0]-2)])
    exec('codes = [Path.MOVETO,'+moves+', Path.CLOSEPOLY]')
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path,facecolor='b',alpha=0.7)
    #plt.text(verts.mean(0)[0],verts.mean(0)[1] , "Cloudband", ha="center", size=14)
    ax = plt.axes()
    ax.add_patch(patch)


def plotblob(m,e,i,datakey,propkey,patch=False):
    '''Helper function to plot contour of blob in plotevent'''
    # Properties of the plot
    plotprops = lineprops[propkey]['blob']
    # Index into mbs array from trkarrs
    if isinstance(e.trkarrs[datakey],dict):
        print 'No track in %s for event %d' %(datakey,e.trkkey)
        return
    if e.trkarrs[datakey].ndim==3:
        for it in xrange((e.trkarrs[datakey].shape[2])):
            ixblob=int(e.trkarrs[datakey][i,1,it])
            if ixblob==-1: continue
            blob = e.blobs[datakey]['mbs'][ixblob]
            contour = e.blobs[datakey]['ch'][ixblob]
            #cX, cY = blob[e.ifld('cX')], blob[e.ifld('cY')]
            cX, cY = e.trkarrs[datakey][i,2,it], e.trkarrs[datakey][i,3,it]
            xverts, yverts = contour[:,0],contour[:,1]
            # Plot centre
            #m.text(cX,cY,str(it),fontsize=12,fontweight='bold')
            #m.plot(cX,cY,'ko',markersize=10)
            exec('m.plot(xverts,yverts,'+plotprops+')')
            if patch:
                addpatch(contour)
    else:
        ixblob=int(e.trkarrs[datakey][i,1])
        if ixblob==-1: return
        blob = e.blobs[datakey]['mbs'][ixblob]
        contour = e.blobs[datakey]['ch'][ixblob]
        cX, cY = e.trkarrs[datakey][i,2], e.trkarrs[datakey][i,3]
        xverts, yverts = contour[:,0],contour[:,1]
        # Plot blob
        #m.text(cX,cY,'dummy',fontsize=12,fontweight='bold')
        #m.plot(cX,cY,'ko',markersize=10)
        exec('m.plot(xverts,yverts,'+plotprops+')')
        if patch:
            addpatch(contour)

def plottrack(m,e,i,datakey,propkey):
    '''Helper function to plot track of blob in plotevent'''
    # Properties of the plot
    plotprops = lineprops[propkey]['trk']
    # Index into mbs array from trkarrs
    if isinstance(e.trkarrs[datakey],dict):
        print 'No track in %s for event %d' %(datakey,e.trkkey)
        return
    trkarr = e.trkarrs[datakey]
    if trkarr.ndim==3:
        for it in xrange((trkarr.shape[2])):
            istrt = np.where(trkarr[:,2,it] != -1)[0][0]
            iend = np.where(trkarr[:,2,it] != -1)[0][-1]
            if i<istrt or i>iend: continue
            tX,tY = trkarr[istrt:i+1,2,it],trkarr[istrt:i+1,3,it]
            cX, cY = trkarr[i,2,it], trkarr[i,3,it]
            # Plot track
            #m.text(cX,cY,str(it),fontsize=12,fontweight='bold')
            m.plot(cX,cY,'ko',markersize=12,markerfacecolor='none',marker='o',linewidth='4')
            exec('m.plot(tX,tY,'+plotprops+')')

    else:
        istrt = np.where(trkarr[:,2] != -1)[0][0]
        iend = np.where(trkarr[:,2] != -1)[0][-1]
        if i<istrt or i>iend: return
        tX,tY = trkarr[istrt:i+1,2],trkarr[istrt:i+1,3]
        cX, cY = trkarr[i,2], trkarr[datakey][i,3]
        # Plot centre
        #m.text(cX,cY,'dummy',fontsize=12,fontweight='bold')
        m.plot(cX,cY,'kx',markersize=12)
        exec('m.plot(tY,tX,'+plotprops+')')


def plotolr(m,lat,lon,olr,compmask=False,compall=False,cbx=0.05):
    '''Helper function to plot pcolor olr for events'''
    plon=lon-1.25;plat=lat+1.25
    x, y = np.meshgrid(plon,plat)
    f = plt.gcf()
    ax = plt.gca()
    if ~np.any(~olr.mask):
        #print 'Nothing to plot in OLR'
        return
    else:
        #olr=olr.data
        #datafig = m.transform_scalar(olr[::-1,:],lon,lat[::-1],len(lon),len(lat))
        #m.imshow(datafig,cmap=plt.cm.PuBu_r,interpolation='nearest')
        #m.pcolor(x,y,olr[:,:],cmap=plt.cm.PuBu_r)
        #m.pcolor(x,y,olr[:,:],cmap=plt.cm.Blues_r)
        m.pcolor(x,y,olr[:,:],cmap=plt.cm.gray)
        if compmask:
            plt.clim(180,215)
            bounds = np.arange(180,215,5)
        if compall:
            plt.clim(210,245)
            bounds = np.arange(215,250,5)
        else:
            plt.clim(100,240)
            bounds = np.arange(140,260,20)
        #vals=np.arange(0,35,2)
        axcol = f.add_axes([cbx,0.3,0.015,0.3])
        plt.colorbar(cax=axcol,boundaries=bounds)
        axcol.yaxis.set_ticks_position('left')
        axcol.yaxis.set_label_position('left')
        plt.axes(ax)

def plotuv(m,lat,lon,u,v,clr,scale=1000.0,lw=.5,keylength=50):
    '''Helper function to plot uv vectors for events'''
    if ~np.any(~u.mask):
        #print 'Nothing to plot in QU, QU'
        return
    else:
        if isinstance(m,tymp): urot,vrot,x,y = m.rotate_vector(u,v,lon,lat,returnxy=True)
        elif isinstance(m,typlt): urot,vrot,x,y = u,v,lon,lat
        mgn = np.ma.sqrt(u**2 + v**2)
        Q = m.quiver(x,y,u,v,color=clr,scale=scale,linewidth=lw,edgecolor=clr)
        plt.quiverkey(Q,0.1,-0.0795,keylength,str(keylength)+' ms$^{-1}$',color=clr,linewidth=lw,labelpos='S')

def plotquv(m,lat,lon,qu,qv,clr,scale=3.0,lw=.5,keylength=0.1):
    '''Helper function to plot qu,qv vectors for events'''
    if ~np.any(~qu.mask):
        #print 'Nothing to plot in QU, QU'
        return
    else:
        if isinstance(m,tymp): urot,vrot,x,y = m.rotate_vector(qu,qv,lon,lat,returnxy=True)
        elif isinstance(m,typlt): urot,vrot,x,y = qu,qv,lon,lat
        mgn = np.ma.sqrt(qu**2 + qv**2)
        Q = m.quiver(x,y,urot,vrot,color=clr,scale=scale,linewidth=lw,edgecolor=clr)
        plt.quiverkey(Q,0.86,-0.08,keylength,str(keylength)+' gkg$^{-1}$s$^{-1}$',color=clr,linewidth=lw,labelpos='S')

def plothgt(m,lat,lon,hgt,lev,colorm=plt.cm.OrRd):
    '''Helper function to plot gph field for events'''
    f = plt.gcf()
    ax = plt.gca()
    if ~np.any(~hgt.mask):
        #print 'Nothing to plot in HGT'
        return
    else:
        plon=lon-1.25;plat=lat+1.25
        x,y = np.meshgrid(plon,plat)
        m.pcolor(x,y,hgt,cmap=colorm)
        if lev==850:
            plt.clim(1480,1520)
            bounds = np.arange(1490,1520,5)
            axcol = f.add_axes([0.91,0.52,0.02,0.25])
            plt.colorbar(cax=axcol,boundaries=bounds)
        if lev==250:
            plt.clim(9900,11000)
            bounds2 = np.arange(10200,10700,100)
            axcol2 = f.add_axes([0.91,0.22,0.02,0.25])
            plt.colorbar(cax=axcol2,boundaries=bounds2)
        plt.axes(ax)

def plotrainstations(m,e,i,rain,lat,lon,datatime,datakey='noaa-olr-0-all',tres=24,threshold=2):
    '''Helper function to plot station rainfall for events'''
    t = e.trkarrstime[datakey][i]
    if tres==24: t = t - (t % 24)
    ix = my.ixtwindow([t],datatime,0)
    if len(ix)==0:print 'No station rainfall available';return
    else: ix=ix[0]
    mask = e.rainmasks[datakey][i,:]
    wetmask = rain[:,ix] > threshold
    inmask = mask & wetmask
    outmask = ~mask & wetmask
    f = plt.gcf()
    ax = plt.gca()
    # Plot data for both in and outside the contour
    if ~np.any(inmask):
        #print 'No rain stations in contour to plot'
        return
    else:
        m.scatter(lon[inmask],lat[inmask],3,rain[inmask,ix],cmap=plt.cm.jet,edgecolor='none')
        plt.clim(threshold,50)
        bounds = np.arange(0,60,10)
        #vals=np.arange(0,35,2)
        axcol = f.add_axes([0.91,0.52,0.02,0.25])
        plt.colorbar(cax=axcol,boundaries=bounds)
        plt.axes(ax)
    if ~np.any(outmask):
        #print 'No rain stations outside contour to plot either'
        return
    else:
        m.scatter(lon[outmask],lat[outmask],3,rain[outmask,ix],cmap=plt.cm.cool,edgecolor='none')
        plt.clim(threshold,50)
        bounds2 = np.arange(0,60,10)
        #vals=np.arange(0,35,2)
        axcol2 = f.add_axes([0.91,0.22,0.02,0.25])
        plt.colorbar(cax=axcol2,boundaries=bounds2)
        plt.axes(ax)

def multiplot(e,m,etime,blobkeylist,trkkeylist,*args):
    '''A plotting multitool that can plot all variables used and produced
    by MetBot package
    USAGE: etime can either be a tuple = (hrsbefore, hrsafter)
                 or a list = [flagday1,flagday2]
           *args are gridded dataset and rainfall station dataset lists
                 [mbskey,data,lat,lon,time,alt_mbskey]
            IMPORTANT: if u,v vector pairs are supplied,
                       they must be supplied in that order
                       if u,v vector pairs for multiple levels are supplied
                       they need to be kept paired'''
    ### UNPACK *ARGS LISTS OF VARIABLES INTO USABLE FORMAT
    hgtlist,uvlist,quvlist,olrlist = [],[],[],[]
    print 'Event:',e.trkkey
    for arg in args:
        v = arg[0];dset, varstr, lev, drv = v.split('-')
        varlev=varstr+'_'+lev
        #print varlev
        exec(varlev+'key,'+varlev+','+varlev+'lat,'+varlev+'lon,'+varlev+'time,'+varlev+'_altkey=arg')
        if varstr=='qv': quvlist.append(['qu'+'_'+lev,varlev])
        elif varstr=='hgt': hgtlist.append(varlev)
        elif varstr=='v': uvlist.append(['u'+'_'+lev,varlev])

    ### PLOT TRKARRS DATA USING INDEX LIST
    masktimes = e.trkarrstime[e.refkey]
    # this because had some issues with some empty trakarrs which won't have the times6hr needed
    for nmbs in e.mbskeys[3:]:
        if ~isinstance(e.trkarrs[nmbs],dict): mbs6hrkey=nmbs
    times6hr = e.trkarrstime[nmbs] # Used for labeling of figure
    dlabels, ixlist, flarr = daylabel(e,etime)
    for i in ixlist:
        # OLR PLOT
        if 'olr_0' in locals():
            eolr = maskeddata(e,i,olr_0,olr_0time,olr_0key)
            plotolr(m,olr_0lat,olr_0lon,eolr)
        # HGT PLOT
        for ph in hgtlist:
            exec('hgt, hgtlat, hgtlon, hgttime, hgtkey = '+ph+', '+ph+'lat, '+ph+'lon, '+ph+'time, '+ph+'key')
            ehgt = maskeddata(e,i,hgt,hgttime,hgtkey)
            plothgt(m,hgtlat,hgtlon,ehgt)
        # STATION RAINFALL PLOT
        if 'rain_0' in locals(): plotrainstations(m,e,i,rain_0,rain_0lat,rain_0lon,rain_0time)
        # BLOB CONTOURS
        for blbkey in blobkeylist: plotblob(m,e,i,blbkey,blobkeylist.index(blbkey))
        # BLOB TRACKS
        for trkkey in trkkeylist: plottrack(m,e,i,trkkey,trkkeylist.index(trkkey))
        # QUV PLOT
        cnt=0
        for pqu, pqv in quvlist:
            exec('qu, qv, qlat, qlon, qtime, qu_altkey, qv_altkey = '+pqu+', '+pqv+', '+pqu+'lat, '+pqu+'lon, '+pqu+'time, '+pqu+'_altkey, '+pqv+'_altkey')
            equ = maskeddata(e,i,qu,qtime,qu_altkey);eqv = maskeddata(e,i,qv,qtime,qv_altkey)
            plotquv(m,qlat,qlon,equ,eqv,cols[cnt]);cnt+=1
        cnt=0
        for pu, pv in uvlist:
            exec('u, v, qlat, qlon, qtime, u_altkey, v_altkey = '+pu+', '+pv+', '+pu+'lat, '+pu+'lon, '+pu+'time, '+pu+'_altkey, '+pv+'_altkey')
            eu = maskeddata(e,i,u,qtime,u_altkey);ev = maskeddata(e,i,v,qtime,v_altkey)
            plotuv(m,qlat,qlon,eu,ev,gcols[cnt]);cnt+=1
        # DRAW EVERYTHING
        redrawmap(m,lns=True)
        # a little thing to handle time labeling relative to event flagging
        tm6,tm24 = times6hr[i],masktimes[i]
        tdate = my.hrnum2dates([tm6])[0]
        eventcode = dlabels[i]
        humandate='Date: %02d-%02d-%02d %02dh00' %(tdate[2],tdate[1],tdate[0],tdate[3])
        plt.title(str(e.trkkey)+': '+eventcode)
        plt.text(30,-49,humandate)
        plt.draw()
        # SAVE FIGURE
        dstring=str(e.trkkey)+'_%02d-%02d-%02d_%02dh00' %(tdate[0],tdate[1],tdate[2],tdate[3])
        fname=figdir+dstring+'.png'
        plt.savefig(fname,dpi=150)
        #sleep(1)
        #raw_input()
        plt.clf()

### MORE CLIMATOLOGICAL TYPE FUNCTIONS
def flagdetailsub(flagarr,ekeys,flagdetail):
    daynos, hr, keys = flagdetail
    icomp = []
    idn, ihr, iky = [], [], []
    ntests = 0
    if daynos:
        ntests+=1
        for d in daynos: ixdn = np.where(flagarr[:,0]==d)[0];idn.extend(list(ixdn))
    if hr: ntests+=1;ihr = np.where(flagarr[:,1]==hr)[0];ihr = list(ihr)
    if keys:
        ntests+=1
        for ik in keys: ixk = np.where(ekeys==ik)[0];iky.extend(list(ixk))
    idn.extend(ihr);idn.extend(iky);ixlist = np.unique(np.asarray(idn))

    ixsub=[]
    if len(ixlist)>0:
        for ix in ixlist:
            if idn.count(ix) >= ntests: ixsub.append(ix)
        ixsub=np.asarray(ixsub)

    return ixsub


def multimasks(s,eventkeys,etime,*args):
    '''Returns a data set with the data and masks for multiple options'''
    # UNPACK THE DATA TUPLES HELD IN ARGS
    varlist = []
    for arg in args:
        v = arg[0];dset, varstr, lev, drv = v.split('-')
        varlev=varstr+'_'+lev
        exec(varlev+'key,'+varlev+','+varlev+'lat,'+varlev+'lon,'+varlev+'time,'+varlev+'_altkey=arg')
        varlist.append(varlev)

    varmasks={}
    for varlev in varlist:
        exec('key, data, datatime, alt_key = '+varlev+'key,'+varlev+','+varlev+'time,'+varlev+'_altkey')
        if isinstance(alt_key,str): key = alt_key
        emdata = np.ma.zeros((0,data.shape[1],data.shape[2]),dtype=np.float32)
        msk = np.ndarray((0,data.shape[1],data.shape[2]),dtype=bool)
        flagdetails = np.ndarray((0,2),dtype=np.float32)
        emdatatime = [];emeventkey = []
        if key=='noaa-olr-0-all':emcXcY = np.ndarray((0,2),dtype=np.float32)

        print "varmasks['%s']: Adding '%s' masks to %s data from %d events" %(varlev,key,varlev,len(eventkeys))
        for k in eventkeys:
            e = s.events[k]
            #print 'Event', k
            # CREATE LIST OF SELECTED INDICES IN TRKARRS AND MASKS
            masktimes = e.trkarrstime[e.refkey]
            dlabels, ixlist, flarr = daylabel(e,etime)
            flagdetails = np.append(flagdetails,flarr,axis=0)
            # COMPUTE MASKS FOR EACH EVENT
            for i in ixlist:
                mdata = maskeddata(e,i,data,datatime,key)
                if np.any(~mdata.mask):mdata.mask[:]=False
                emdata = np.append(emdata,mdata.data[np.newaxis,:,:],axis=0)
                msk     = np.append(msk , mdata.mask[np.newaxis,:],axis=0)
                if isinstance(e.trkarrs[key],np.ndarray):emdatatime.append(e.trkarrstime[key][i])
                elif isinstance(e.trkarrs[key],dict):emdatatime.append(masktimes[i])
                emeventkey.append(k)

                if key=='noaa-olr-0-all':emcXcY = np.append(emcXcY,e.trkarrs[key][i,2:][np.newaxis,:],axis=0)
        emdata.mask = msk
        varmasks[varlev] = (emdata,np.asarray(emdatatime),np.asarray(emeventkey),flagdetails)
        if key=='noaa-olr-0-all':varmasks[varlev] = (emdata,np.asarray(emdatatime),np.asarray(emeventkey),flagdetails,emcXcY)

    return varmasks

def multirainmasks(s,eventkeys,etime,*args,**kwargs):
    '''Returns a data set with the data and masks for multiple options'''
    # UNPACK THE DATA TUPLES HELD IN ARGS
    varlist = []
    for arg in args:
        v = arg[0];dset, varstr, lev, drv = v.split('-')
        exec(varstr+'key,'+varstr+','+varstr+'lat,'+varstr+'lon,'+varstr+'time,'+varstr+'_altkey=arg')
        varlist.append(varstr)
    dailyonly=False
    for key in kwargs:
        dailyonly=kwargs[key]
    rainmasks={}
    for varstr in varlist:
        exec('key, data, datatime, alt_key = '+varstr+'key,'+varstr+','+varstr+'time,'+varstr+'_altkey')
        #print varstr,data.shape
        ermdata = np.ma.zeros((0,data.shape[0]),dtype=np.float32)
        msk = np.ndarray((0,data.shape[0]),dtype=bool)
        flagdetails = np.ndarray((0,2),dtype=np.float32)
        ermdatatime = [];ermeventkey=[]
        print "rainmasks['%s']: Adding '%s' masks to %s data from %d events" %(key,key,varstr,len(eventkeys))
        for k in eventkeys:
            e = s.events[k]
            #print 'Event', k
            # CREATE LIST OF SELECTED INDICES IN TRKARRS AND MASKS
            masktimes = e.trkarrstime[e.refkey]
            dlabels, ixlist, flarr = daylabel(e,etime)
            if dailyonly:
                ix24 = np.where(flarr[:,1]%24==0)[0]
                dlabels,ixlist,flarr=dlabels[ix24],ixlist[ix24],flarr[ix24,:]
            flagdetails = np.append(flagdetails,flarr,axis=0)
            for i in ixlist:
                rmdata = rainmaskeddata(e,i,data,datatime,alt_key)
                if isinstance(rmdata,bool):
                    print 'Finished available dates for',key,'data set'
                    continue
                ermdata = np.append(ermdata,rmdata.data[np.newaxis,:],axis=0)
                msk     = np.append(msk ,   rmdata.mask[np.newaxis,:],axis=0)
                ermdatatime.append(e.trkarrstime[alt_key][i])
                ermeventkey.append(k)

        ermdata.mask = ~msk
        ermdatatime, ermeventkey = np.asarray(ermdatatime),np.asarray(ermeventkey)
        rainmasks[varstr] = (ermdata,ermdatatime,ermeventkey,flagdetails)
        if False:
            pickf = open('/home/neil/work/computervision/metblobs/SA/pickles/rainmasks.msk','w')
            blb.cPickle.dump(rainmasks,pickf)
            pickf.close()

    return rainmasks

def composite(varmask,kv,flagdetail,freqthresh=False):
    '''This function will return the composited data
    USAGE: varmasks (dict) returned by multimasks
           kv is the key into the varmask you want
           flagdetail (tuple) of (daynos,hr,key)
                                       dayno is the event day number of events
                                       hr is hour of event or before event
                                       key can be a list of keys of events'''
    print 'Calculating composite event for',kv
    if not kv=='olr_0':
        exec('data,datatime,ekey,flagarr = varmask')
    elif kv=='olr_0':
        exec('data,datatime,ekey,flagarr,cxcy = varmask')
    ixsub = flagdetailsub(flagarr,ekey,flagdetail)
    if len(ixsub)>0:
        data,datatime,ekey,flagarr = data[ixsub,:,:],datatime[ixsub],ekey[ixsub],flagarr[ixsub,:]
        if kv=='olr_0': cxcy=cxcy[ixsub,:]
    datamn = np.where(data.mask,np.nan,data.data);datamn = sps.nanmean(datamn,0)
    #datamn = sps.nanmean(data.data,0)
    if freqthresh:
        nfreq = np.where(data.mask,0.0,1.0)
        nfreq = np.sum(nfreq,0)/nfreq.shape[0]
        nfreqmsk = ~ (nfreq >= freqthresh)
        out = (np.ma.MaskedArray(datamn,mask=nfreqmsk,dtype=np.float32), nfreq, ixsub)
    else:
        out = (np.ma.MaskedArray(datamn,dtype=np.float32), 0, ixsub)

    return out

def relativegrid(cxcy,lon,lat,data,ycentre=False):
    '''Regrid data on a grid centred on xycentre
    Generally each xycentre corresponds to one timestep
    So usage can be: xycentre (time x 2, ndarray)
                     data      (time x lat x lon)
                     ycentre = False by default, so only create data relative to longitudinal centres

    NOTE: This will only work on a regularly spaced lat x lon grid'''

    # TEST THAT DATA IS LAT X LON
    nln, nlt = len(lon), len(lat)
    nt,ny,nx = data.shape
    if ny!=nlt or nx!=nln: print "Sure your data is in time x lat x lon? I don't thinks so!" ;return

    # CREATE RELATIVE LATXLON COORDS AND THE DATA GRID THEY DESCRIBE
    dlt=(lat[1:]-lat[:-1])[0];rlat=lat[-1]-lat[0];rellat = np.arange((rlat*-1),rlat+dlt,dlt)
    dln=(lon[1:]-lon[:-1])[0];rlon=lon[-1]-lon[0];rellon = np.arange((rlon*-1),rlon+dln,dln)
    rln = len(rellon)
    if ~ ycentre: rlt = len(lat)
    else: rlt = len(rellat)
    relgrid = np.zeros((nt,rlt,rln),dtype=np.float32);relgrid[:]=np.nan
    relgridmsk = np.ones((nt,rlt,rln),dtype=np.bool)

    # MAIN LOOP
    for i in xrange(nt):
        cx,cy = cxcy[i]; rdata = data[i,:,:]
        if cx==-1.0 and cy==-1: print 'No data in timeslice=',i;continue
        xoff, yoff = lon-cx,lat-cy
        ix = np.where(rellon>xoff[0])[0][0]
        iy = np.where(rellat<yoff[0])[0][0]
        if ~ ycentre:
            relgrid[i,:,ix:ix+nx] = rdata.data
            relgridmsk[i,:,ix:ix+nx] = rdata.mask
        else:
            relgrid[i,iy:iy+ny,ix:ix+nx] = rdata.data
            relgridmsk[i,iy:iy+ny,ix:ix+nx] = rdata.mask

    relgrid = np.ma.MaskedArray(relgrid,mask=relgridmsk,dtype=np.float32)

    # FINALLY WE CUT THE GRID SMALLER AGAIN
    if ~ ycentre:rellat=lat; iy1,iy2 = 0,len(rellat)
    else: iy1,iy2 = int(rlt*0.25),int(rlt*0.75)+1

    ix1,ix2 = int(rln*0.25),int(rln*0.75)+1

    return relgrid[:,iy1:iy2,ix1:ix2], rellat[iy1:iy2], rellon[ix1:ix2]

def eventsremovedcomp(varmasks,kv,rawdata,rawdatatime,flagdetail):
    '''similar to composite except meant to remove influence of events from raw data'''

    if not kv=='olr_0': exec('data,datatime,ekey,flagarr = varmasks["'+kv+'"]')
    elif kv=='olr_0': exec('data,datatime,ekey,flagarr,cxcy = varmasks["'+kv+'"]')

    ixsub = flagdetailsub(flagarr,ekey,flagdetail)
    if len(ixsub)>0:
        data,datatime,ekey,flagarr = data[ixsub,:,:],datatime[ixsub],ekey[ixsub],flagarr[ixsub,:]
        if kv=='olr_0': cxcy=cxcy[ixsub,:]

    rawfreq = np.where(rawdata.mask,0.0,1.0)
    rawfreq = np.sum(rawfreq,0)/rawfreq.shape[0]

    rawmaskrem = rawdata.mask.copy()
    ixdatasub=[]
    for hr in datatime:
        ixt = np.where(rawdatatime==hr)[0]
        if len(ixt)>0: ixdatasub.append(ixt[0])
    ixdatasub = np.unique(ixdatasub)
    rawmaskrem[ixt,:,:] = True
    rawfreqrem = np.where(rawmaskrem,0.0,1.0)
    rawfreqrem = np.sum(rawfreqrem,0)/rawfreqrem.shape[0]

    #nfreq = np.where(data.mask,0.0,1.0)
    #nfreq = np.sum(nfreq,0)/nfreq.shape[0]
    #nfreqmsk = ~ (nfreq >= freqthresh)

    datamn = np.where(rawdata.mask,np.nan,rawdata.data)
    datamn = sps.nanmean(datamn,0)

    datamnrem = np.where(rawmaskrem,np.nan,rawdata.data)
    datamnrem = sps.nanmean(datamnrem,0)

    return datamn, datamnrem, np.asarray(ixdatasub)

def speedup(s,kchoice,raintup,rlist,idx):
    '''This is the function run by threadedmultimasks'''
    print 'Running thread...'
    rainmask = multirainmasks(s,kchoice[idx[0]:idx[1]],(-48,48),raintup,dailyonly=True)
    rlist.append(rainmask)

def threadedmultirainmasks(s,kchoice,raintup,tms=(-48,48)):
    '''simple threaded version of multirainmask to speed up things
    couple of things to be completed to implement nthreads=(any other no.)'''
    # Get the threads running
    nthreads=4
    nl=len(kchoice);nt=nl/nthreads;ixx=[0,nt]
    for ith in range(0,nthreads):
        idx=[ixx[0]+nt*ith,ixx[1]+nt*ith]
        if ith==(nthreads-1):idx=[ixx[0]+nt*ith,nl]
        exec('ls'+str(ith)+'=[]')
        exec('thread'+str(ith)+' = threading.Thread(target=speedup,args=(s,kchoice,raintup,ls'+str(ith)+',idx))')
        exec('thread'+str(ith)+'.start()')
    # dirty way to wait for all threads to finish
    while len(ls3)==0 or len(ls2)==0 or len(ls1)==0 or len(ls0)==0:
        dum="why are we waiting, going to keep on waiting"
    # pull it all together
    if len(ls3)>0 and len(ls2)>0 and len(ls1)>0 and len(ls0)>0:
        print 'Success, now concatentating'
        maskdata, maskhrtime, maskkeys, maskflagarr= ls0[0]['rain']
        maskd, maskmsk = maskdata.data, maskdata.mask
        for i in [1,2,3]:
            exec('moredata, morehrtime, morekeys, moreflagarr = ls'+str(i)+'[0]["rain"]')
            maskd = np.append(maskd,moredata.data,axis=0)
            maskmsk = np.append(maskmsk,moredata.mask,axis=0)
            maskhrtime = np.append(maskhrtime,morehrtime,axis=0)
            maskkeys = np.append(maskkeys,morekeys,axis=0)
            maskflagarr = np.append(maskflagarr, moreflagarr,axis=0)
        maskdata = np.ma.MaskedArray(maskd,mask=maskmsk)
        rainmasks={'rain':(maskdata,maskhrtime,maskkeys,maskflagarr)}
    else:
        print "Oops, something didn't work here"

    return rainmasks
