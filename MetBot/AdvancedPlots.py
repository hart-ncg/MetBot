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
    if mbkl[0]=='um':dclim=(1,9,1)
    if mbkl[0]=='cmip5':dclim=(1,9,1)

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season=='dryseason':mns=[4,5,6,7,8,9]
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
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
        my.xtickfonts();my.ytickfonts()
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(file_suffix+'_flagonly.png',dpi=150)
        else: plt.savefig(file_suffix+'.png',dpi=150)

    return msklist

def spatiofreq3_season(s,lat,lon,yrs,eventkeys,figno=1,season='coreseason',\
     key='noaa-olr-0-all',res='noaa',dclim=(1,7,1),flagonly=False,file_suffix='test',savefig=False,\
     fontdict=False):
    '''spatiofreq3_season(s,lat,lon,yrs,figno=1,season='coreseason',\
                          flagonly=False,file_suffix='test')
    Produces subplots of cloud-band gridpoint count by month
    spatiofreq3 similar to spatiofreq2 but with option by RJ to fix res to noaa
    '''
    if not fontdict:fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbkl=key.split('-')

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season=='dryseason':mns=[4,5,6,7,8,9]
    elif isinstance(season,list):
        mns=season

    if mbkl[0]=='noaa':
        m, f = blb.SAfrBasemap(lat[6:-8],lon[3:-2],drawstuff=True,prj='cyl',fno=figno,rsltn='l')
    else:
        if res=='native':
            m, f = blb.SAfrBasemap(lat[13:-27],lon[6:-5],drawstuff=True,prj='cyl',fno=figno,rsltn='l')
        else:
            m, f = blb.SAfrBasemap(lat[6:-8], lon[3:-2], drawstuff=True, prj='cyl', fno=figno,
                                   rsltn='l')  # lat and lon specified by number - diff for UM at native res

    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])


    if len(mns)==12:plt.close();plt.figure(figsize=[12,10])
    elif len(mns)==6:plt.close();plt.figure(figsize=[13,8])
    cnt=1
    msklist=[]
    for mn in mns:
        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)
        if flagonly:
            allmask=stats.spatiofreq3(m,s,lat,lon,yrs,eventkeys,\
                    clim=dclim,month=mn,flagonly=True,fontdict=fd)
        else:
            allmask=stats.spatiofreq3(m,s,lat,lon,yrs,eventkeys,clim=dclim,\
                                      month=mn,flagonly=False,fontdict=fd)
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
        my.xtickfonts();my.ytickfonts()
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(file_suffix+'_res-'+res+'_flagonly_spatiofreq.png',dpi=150)
        else: plt.savefig(file_suffix+'_res-'+res+'_spatiofreq.png',dpi=150)

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


def gridrainmap_season(s,raingrid,rain,lat,lon,rdtime,eventkeys,yrs,figno=1,season='coreseason',key='noaa-olr-0-0',ptype='per_ttt',file_suffix='test',savefig=False):
    '''Produces subplots of ttt rainfall by month
    need to open the rain data with lon and lat and also the synop file
    see e.g. rainmap_test_months.py

    plot types
        tot_all - sum of total precip
        tot_ttt - sum of precip under TTTs
        per_ttt - percent of precip under TTTs
    '''
    mbkl=key.split('-')
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    if yrs == None:
        print "No years specified"
        yrs = np.unique(rdtime[:,0])
    else:
        print "years specified"
        print yrs
        yrs=yrs
    print yrs

    nlon=len(lon)
    nlat=len(lat)

   # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season == 'dryseason':mns = [4, 5, 6, 7, 8, 9]
    elif isinstance(season,list):
        mns=season

    # Draw basemap
    m, f = blb.SAfrBasemap(lat,lon,drawstuff=True,prj='cyl',fno=figno,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    cnt=1
    msklist=[]
    for mn in mns:
        print mn

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data (allmask)
        speckeys = stats.specificmon(s,eventkeys,yrs,mn,mbkl[0])

        #All rain
        firstyear=yrs[0]
        lastyear=yrs[len(yrs)-1]
        raindat=np.where((rdtime[:,0]>=firstyear) & (rdtime[:,0]<=lastyear))
        rainperiod=rain[raindat,:,:]
        rainperiod=np.squeeze(rainperiod)
        newdates=rdtime[raindat]
        ix=np.where((newdates[:,1]==mn))
        rainmon=rainperiod[ix,:,:]
        rainmon=np.squeeze(rainmon)
        rainsum_all=np.nansum(rainmon,0)

        if ptype=='tot_all':
            allmask=rainsum_all
        else:
            #TTT rain
            rain4plot = stats.griddedrainmasks(s,speckeys,raingrid,refkey=key)

            nevent=rain4plot.shape[0]
            data2sum=np.zeros((nevent,nlat,nlon),dtype=np.float32)
            for e in xrange(nevent):
                data=rain4plot[e]
                ntstep=len(data[:,0,0])
                data2sum[e,:,:]=np.nansum(data,0)

            rainsum_ttt=np.nansum(data2sum,0)

            if ptype=='tot_ttt':
                allmask=rainsum_ttt
            elif ptype=='per_ttt':
                rainperc_ttt=(rainsum_ttt/rainsum_all)*100
                allmask=rainperc_ttt
            elif ptype=='rainperttt':
                dclim=(1,9,1)
                spatiofreq=stats.spatiofreq2(m,s,lat,lon,yrs,speckeys,clim=dclim,key=key,month=mn,flagonly=True)

        #Plot
        plon,plat = np.meshgrid(lon,lat)
        if ptype=='per_ttt':
            clevs=[0,10,20,30,40,50,60,70,80]
            #cm=plt.cm.viridis
            cm = plt.cm.YlGnBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        elif ptype=='tot_all':
            clevs=[0,800,1600,2400,3200,4000,4800,5600]
            cm=plt.cm.viridis
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        elif ptype=='tot_ttt':
            #clevs=[0,150,300,450,600,750,900,1050,1200,1350,1500]
            clevs=[0,250,500,750,1000,1250,1500,1750,2000]
            cticks = [0,500,1000,1500,2000]
            #cm=plt.cm.viridis
            #cm=plt.cm.Blues
            #cm=plt.cm.terrain_r
            #cm=plt.cm.GnBu
            cm=plt.cm.YlGnBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        plt.title(monthname[mn-1])

        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.02,wspace=0.02,hspace=0.2)
    #axcl=g.add_axes([0.9, 0.15, 0.02, 0.7])
    axcl = g.add_axes([0.91, 0.1, 0.01, 0.12]) # small cbar for paper
    cbar = plt.colorbar(cs,cax=axcl,ticks=cticks)
    if ptype=='per_ttt':cbar.set_label('%')
    else:cbar.set_label('mm')

    if savefig:
        plt.savefig('Grid_pcent_TTTrain_months_'+file_suffix+'.png',dpi=150)

    return msklist

def gridrainmap_bias_season(s,raingrid,rain,lat,lon,rdtime,eventkeys,yrs,figno=1,season='coreseason',key='um-olr-0-0',ptype='diff_tot',file_suffix='test',savefig=False):
    '''Produces subplots of ttt rainfall by month
    need to open the rain data with lon and lat and also the synop file
    see e.g. rainmap_test_months.py

    plot types
        diff_tot - UM total rain - TRMM total rain
        diff_ttt - UM ttt rain - TRMM ttt rain
        percent - (diff_ttt/diff_tot)*100 - % of bias from ttt bias

    '''
    mbkl=key.split('-')
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    if yrs == None:
        print "No years specified"
        yrs = np.unique(rdtime[:,0])
    else:
        print "years specified"
        print yrs
        yrs=yrs
    print yrs

    nlon=len(lon)
    nlat=len(lat)

    # Open TRMM data
    refdir='/ouce-home/students/bras2392/work/impala/data/refdata/'
    trmmfile=refdir+'trmm_3b42v7_n216.daily.precip.1998_2013.mon.nc'

    train, time, tlat, tlon, tdtime = \
        mync.opentrmm(trmmfile,'pcp','trmm',subs='SA_TR')
    tdtime[:,3]=0 # editing the hour time in date
    traingrid=(train,tdtime,(tlon,tlat))

    tnlon=len(tlon)
    tnlat=len(tlat)

    # Open satellite events
    satdir='/ouce-home/students/bras2392/PycharmProjects/ttcb_um_test/test/'
    satsyfile = satdir + "NOAA-OLR.synop"
    sats = sy.SynopticEvents((),[satsyfile],COL=False)
    satks = sats.events.keys();satks.sort()
    tkey='noaa-olr-0-0'
    tmbkl=tkey.split('-')

   # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season

    # Draw basemap
    m, f = blb.SAfrBasemap(lat,lon,drawstuff=True,prj='cyl',fno=figno,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    cnt=1
    msklist=[]
    for mn in mns:
        print mn

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data (allmask)
        speckeys = stats.specificmon(s,eventkeys,yrs,mn,mbkl[0])
        tspeckeys = stats.specificmon(sats,satks,yrs,mn,tmbkl[0])

        #All rain
        firstyear=yrs[0]
        lastyear=yrs[len(yrs)-1]
        raindat=np.where((rdtime[:,0]>=firstyear) & (rdtime[:,0]<=lastyear))
        traindat=np.where((tdtime[:,0]>=firstyear) & (tdtime[:,0]<=lastyear))

        rainperiod=rain[raindat,:,:]
        trainperiod=train[traindat,:,:]

        rainperiod=np.squeeze(rainperiod)
        trainperiod=np.squeeze(trainperiod)

        newdates=rdtime[raindat]
        tnewdates=tdtime[traindat]

        ix=np.where((newdates[:,1]==mn))
        tix=np.where((tnewdates[:,1]==mn))

        rainmon=rainperiod[ix,:,:]
        trainmon=trainperiod[tix,:,:]

        rainmon=np.squeeze(rainmon)
        trainmon=np.squeeze(trainmon)

        rainsum_all=np.nansum(rainmon,0)
        trainsum_all=np.nansum(trainmon,0)

        allrainsum_diff=rainsum_all-trainsum_all

        if ptype=='diff_tot':
            allmask=allrainsum_diff
        else:
            #TTT rain
            rain4plot = stats.griddedrainmasks(s,speckeys,raingrid,refkey=key)
            train4plot = stats.griddedrainmasks(sats,tspeckeys,traingrid,refkey=tkey)

            nevent=rain4plot.shape[0]
            tnevent=train4plot.shape[0]

            data2sum=np.zeros((nevent,nlat,nlon),dtype=np.float32)
            tdata2sum=np.zeros((tnevent,tnlat,tnlon),dtype=np.float32)

            for e in xrange(nevent):
                data=rain4plot[e]
                ntstep=len(data[:,0,0])
                data2sum[e,:,:]=np.nansum(data,0)

            for e in xrange(tnevent):
                data=train4plot[e]
                ntstep=len(data[:,0,0])
                tdata2sum[e,:,:]=np.nansum(data,0)

            rainsum_ttt=np.nansum(data2sum,0)
            trainsum_ttt=np.nansum(tdata2sum,0)

            tttrainsum_diff=rainsum_ttt-trainsum_ttt

            if ptype=='diff_ttt':
                allmask=tttrainsum_diff
            elif ptype=='percent':
                biasperc_ttt=(tttrainsum_diff/allrainsum_diff)*100
                allmask=biasperc_ttt
            elif ptype=='rainperttt':
                dclim=(1,9,1)
                spatiofreq=stats.spatiofreq2(m,s,lat,lon,yrs,speckeys,clim=dclim,key=key,month=mn,flagonly=True)

        #Plot
        plon,plat = np.meshgrid(lon,lat)
        if ptype=='diff_tot':
            clevs=[-1500,-1000,-500,0,500,1000,1500]
            cm=plt.cm.RdBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        elif ptype=='diff_ttt':
            clevs=[-450,-300,-150,0,150,300,450]
            cm=plt.cm.RdBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        elif ptype=='percent':
            clevs=[0,10,20,30,40,50,60,70,80,90,100]
            #Mask negative values
            allmask=np.ma.masked_where(allmask < 0, allmask)
            cm=plt.cm.RdPu
            cm.set_bad(color='white')
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='max')
        else:
            cm=plt.cm.RdYlBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,cmap=cm)
            plt.colorbar()
        plt.title(monthname[mn-1])

        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.02,wspace=0.02,hspace=0.2)
    axcl=g.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs,cax=axcl)
    if ptype=='percent':cbar.set_label('%')
    else:cbar.set_label('mm')

    if savefig:
        plt.savefig('Grid_bias_TTTrain_months_'+file_suffix+'.png',dpi=150)

    return msklist