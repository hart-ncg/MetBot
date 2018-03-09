'''EventStats.py: A module of MetBot package

   Main Task is to calculate and plot temporal characteristics of events'''
import SynopticAnatomy as sy
import numpy as np
import matplotlib.pyplot as plt
from mytools import points_inside_poly
import scipy.stats as st
import mynetcdf as mync
import mytools as my
import MetBlobs as blb
import datetime
import cPickle
from time import time as timer
season=[8,9,10,11,12,1,2,3,4,5,6,7]
coreseason=[10,11,12,1,2,3]
monthlist=[1,2,3,4,5,6,7,8,9,10,11,12]
monthends = [31,28,31,30,31,30,31,31,30,31,30,31]
monthends360 = [30,30,30,30,30,30,30,30,30,30,30,30]
monthends_leap = [31,29,31,30,31,30,31,31,30,31,30,31]
monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']
monthstrseason=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']
monthstr_calendar=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec']
mndict={}
for i in xrange(len(season)): mndict[season[i]]=monthstr[i]


domains={'SA': ((-60.0,0.0),(0.0,80.0)),
'Fauch08': ((-40.0,-15.0),(7.5,70.0),),
'SH': ((-90.0,0.0),(0,360)),}
sbst='SA'
vlt, vln = domains[sbst][0], domains[sbst][1]

SPCPyears = np.asarray([2004,2009])
CPyears = np.asarray([1987,1990,1994,2002, 2004, 2009])
EPyears = np.asarray([1982,1991,1997])
ELyears = np.asarray([1982,1986,1987,1991,1994,1997,2002,2004,2009])
LAyears = np.asarray([1984,1988,1995,1998,1999,2000,2007,2010,2011])
nonENSOyrs=range(1979,2013)
for y in np.hstack((LAyears,ELyears)):
    nonENSOyrs.remove(y)
nonENSOyears=np.asarray(nonENSOyrs)
dudxyears=np.asarray([1980,1981,1983,1986,1987,1990,1993,1995,1996,1999,2001,2003,2004,2007,2008])
nodudxyears=np.asarray([2000,1992,1984])
eastdudxyears=np.asarray([1979,1982,1985,1988,1989,1991,1994,1997,1998,2002,2006,2009])

#figdir='/home/neil/work/computervision/metblobs/SA/statsfigs/'
#rainfigdir='/home/neil/work/computervision/metblobs/SA/rain/'

def spatialsubset(s,eventkeys,cutlon=45.0):
    '''Get tracks from subsetted domain
    Currently just splits tracks left and right of given lon'''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    x,y=np.array(()),np.array(())
    westkeys, eastkeys = [], []
    for k in eventkeys:
        e = s.events[k]
        #x,y = e.trkcX.mean(), e.trkcY.mean()
        x,y = e.trkcX[0], e.trkcY[0]
        if x < cutlon:
            westkeys.append(k)
        else:
            eastkeys.append(k)
    return westkeys, eastkeys

def timesubset(s,eventkeys,edates,units,cal):
    '''subsetofkeys = timesubset(s,eventkeys,edates)

    Get tracks from certain months or years
    This is done by testing whether the trkarrstimes of an event match any of 
    the dates specified in edates or if edates is list of dates.
    
    ALTERNATIVELY this is done by checking whether these trkarrastimes fall
    between the dates in edates when specified in a tuple (date1, date2)
    where date1/date2 are each lists of specifying [YYYY,MM,DD,HH]
    '''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    e = s.events[eventkeys[0]]
    refkey=e.refkey
    dlist = refkey.split('-')

    keysubset=[]
    if isinstance(edates,list):
        if 'noaa' in dlist:
            hrs = my.dates2hrnum(edates)
        elif 'hadam3p' in dlist:
            hrs = my.dates2hrnum(edates,\
                    units="hours since 1959-12-01 00:00:00",calendar="360_day")
            hrs=hrs/24
        elif 'um' in dlist:
            hrs = my.dates2hrnum(edates, units="hours since 1978-09-01 00:00:00", calendar="360_day");
            hrs = hrs / 24
        for k in eventkeys:
            e = s.events[k]
            time24hr = e.trkarrstime[refkey]
            for hr in hrs:
                if np.any(hr==time24hr): keysubset.append(k)
    elif isinstance(edates,np.ndarray):
        if 'noaa' in dlist:
            hrs = my.dates2hrnum(edates)
        elif 'hadam3p' in dlist:
            hrs = my.dates2hrnum(edates,\
                    units="hours since 1959-12-01 00:00:00",calendar="360_day")
            hrs=hrs/24
        elif 'um' in dlist:
            hrs = my.dates2hrnum(edates, units="hours since 1978-09-01 00:00:00", calendar="360_day");
            hrs = hrs / 24
        for k in eventkeys:
            e = s.events[k]
            time24hr = e.trkarrstime[refkey]
            for hr in time24hr:
                if np.any(hrs==hr): keysubset.append(k)
    elif isinstance(edates,tuple):
        # Edited this part to be flexible to model - RJ 2017
        # if 'noaa' in dlist:
        #     hrs = my.dates2hrnum(edates)
        #     #hrs=hrs/24 - need this for noaa cdr but not noaa noaa
        # elif 'hadam3p' in dlist:
        #     hrs = my.dates2hrnum(edates,\
        #             units="hours since 1959-12-01 00:00:00",calendar="360_day")
        #     hrs=hrs/24
        # elif 'um' in dlist:
        #     hrs = my.dates2hrnum(edates, \
        #             units="hours since 1978-09-01 00:00:00", calendar="360_day");
        #     hrs = hrs / 24
        hrs = my.dates2hrnum(edates,units,cal)
        #Check the time unit
        k=eventkeys[0]
        e = s.events[k]
        time24hr = e.trkarrstime[refkey]
        num_digits=len(str(int(time24hr[0])))
        print num_digits
        if num_digits==5:
            hrs=hrs
        elif num_digits==7: # actually Im not sure if this is needed now that reading in units
            hrs=hrs/24
        else:
            print 'Unknown time listing'
        hrmn, hrmx = hrs
        for k in eventkeys:
            e = s.events[k]
            time24hr = e.trkarrstime[refkey]
            if np.any((time24hr>=hrmn) & (time24hr<=hrmx)): keysubset.append(k)
    keysubset = list(np.unique(np.asarray(keysubset)))

    return keysubset

def specificseasons(s,eventkeys,seasonstartyr,startd=[10,01],endd=[4,30]):
    '''Uses timesubset function to return keys from specific (August-July)
    seasons specified by their starting year'''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    specifickeys=[]
    for styr in seasonstartyr:
        dstart= [styr,startd[0],startd[1],0]
        dend = [styr+1,endd[0],endd[1],0]
        datetup = (dstart,dend)
        seasevents = timesubset(s,eventkeys,datetup)
        specifickeys.extend(seasevents)

    return specifickeys

def specificmon(s,eventkeys,yrs,month,units,cal):
    '''Uses timesubset function to return keys for specific month and years'''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    specifickeys=[]
    for yr in yrs:
        dstart=[yr,month,1,0]
        if cal=='360_day':
            lastday=monthends360[month-1]
        else:
            lastday=monthends[month-1]
        dend = [yr,month,lastday,0]
        datetup = (dstart,dend)
        seasevents = timesubset(s,eventkeys,datetup,units,cal)
        specifickeys.extend(seasevents)

    return specifickeys

def getwinterdates(s,ks,yrs,dset):
    mons=[6,7,8]
    key=dset+'-olr-0-0'
    datelist=[]
    for m in mons:
        keys=specificmon(s,ks,yrs,m,dset)
        for k in keys:
            e=s.events[k]
            hrs=e.trkarrstime[key]
            dates=my.hrnum2dates(hrs,units="hours since 1800-01-01 00:00:0.0",calendar='gregorian')
            datelist.append(dates)

    return datelist

def seasonalcycle(s,eventkeys,years=False,season=[8,9,10,11,12,1,2,3,4,5,6,7]):
    '''Calculate seasonal cycle of cloudband frequency
    Can specify eventkeys is specified False if want all events
    Can specify years or allow autodiscovery'''
    mbskeys = s.mbskeys
    refkey = mbskeys[0]
    basekey = mbskeys[1]
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    # GET CENTRAL TIME OF EVENT
    edts=[]
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[refkey]['mbt'][e.ixflags]
        if len(dts)>1:
            dts = dts[len(dts)/2]
        else:
            dts=dts[0]
        edts.append(dts)
    edts = np.asarray(edts)
    # BUILD SEASONAL CYCLES OF EACH YEAR
    yrs = np.unique(edts[:-1,0])[:-1];dirtyfix=False
    ###### WHAT IS THIS DIRTYFIX ABOUT??? ##################
    if len(yrs)==0:yrs = np.unique(edts[:-1,0]);dirtyfix=True
    ############## HUH? ##############################
    if isinstance(years,list) or isinstance(years,np.ndarray): yrs=years
    scycle = np.zeros((len(yrs),12))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        for imn in xrange(len(season)):
            mn=season[imn]
            if imn < 5:
                ix = np.where((edts[:,0]==yr) & (edts[:,1]==mn))[0]
            elif dirtyfix:
                ix = np.where((edts[:,0]==yr) & (edts[:,1]==mn))[0]
            else:
                ix = np.where((edts[:,0]==yr+1) & (edts[:,1]==mn))[0]
            scycle[iyr,imn] = len(ix)
    cyclestats=np.hstack((scycle.mean(0)[:,np.newaxis],\
                          scycle.std(0)[:,np.newaxis]))
    return scycle, cyclestats, yrs

def scycle_rainfall(s,eventkeys,raindset='wrc',years=False,\
                    season=[8,9,10,11,12,1,2,3,4,5,6,7]):
    '''Calculate seasonal cycle of mean cloudband intensity
    Can specify eventkeys is specified False if want all events
    Can specify years or allow autodiscovery'''
    mbskeys = s.mbskeys
    refkey = mbskeys[0]
    basekey = mbskeys[1]
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    # GET CENTRAL TIME OF EVENT, MEAN RAINFALL AND TOTAL RAINFALL
    # RJ - note total rainfall is based on mean for each blob in the event
    edts=[]
    emrain=np.ndarray((0,8),dtype=np.float32)
    etrain=np.ndarray((0,8),dtype=np.float32)
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[refkey]['mbt'][e.ixflags]
        #dts = s.blobs[basekey]['mbt'][e.trk]
        if len(dts)>1:
            dts = dts[len(dts)/2]
        else:
            dts=dts[0]

        if len(e.rainfall[raindset])==0:
            print 'No data availble in',raindset,'for event',k
            continue

        edts.append(dts)
        mrain = st.nanmean(e.rainfall[raindset],0)
        mrain = np.append(mrain,e.rainfall[raindset].shape[0])
        emrain = np.append(emrain,mrain[np.newaxis,:],axis=0)

        train = np.nansum(e.rainfall[raindset],0)
        if isinstance(train,np.float):
            print 'bust'
            train=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
        train = np.append(train,e.rainfall[raindset].shape[0])
        etrain = np.append(etrain,train[np.newaxis,:],axis=0)
    edts = np.asarray(edts)


    # BUILD SEASONAL CYCLES FROM ALL YEARS
    yrs = np.unique(edts[:-1,0])[:-1]
    if years: yrs=years
    scycle_mrain = np.zeros((len(yrs),12,8))
    scycle_train = np.zeros((len(yrs),12,8))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        for imn in xrange(len(season)):
            mn=season[imn]
            if imn < 5:
                ix = np.where((edts[:,0]==yr) & (edts[:,1]==mn))[0]
            else:
                ix = np.where((edts[:,0]==yr+1) & (edts[:,1]==mn))[0]
            scycle_mrain[iyr,imn,:] = st.nanmean(emrain[ix,:])
            scycle_train[iyr,imn,:] = np.nansum(etrain[ix,:])


    return scycle_mrain, scycle_train, yrs


def scycle_rainsum(s,eventkeys,raindset='wrc',years=False,heavy=False,season=[8,9,10,11,12,1,2,3,4,5,6,7]):
    '''Calculate seasonal cycle of rainfall contributed by TTCBs
    The total rainfall under an OLR flag summed
    RJ July 2016 based on scycle_rainfall
    Can specify eventkeys is specified False if want all events
    heavy key word to sum only over heavy days (thres prescribed in addeventrain)
    Can specify years or allow autodiscovery'''
    if heavy:
        rainind=8
    else:
        rainind=7
    mbskeys = s.mbskeys
    refkey = mbskeys[0]
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    # GET CENTRAL TIME OF EVENT, TOTAL RAINFALL
    edts=[]
    etrain=[]
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[refkey]['mbt'][e.ixflags]
        if len(dts)>1:
            dts = dts[len(dts)/2] # Picks the middle date
        else:
            dts=dts[0]

        if len(e.rainfall[raindset])==0:
            print 'No data availble in',raindset,'for event',k
            continue

        edts.append(dts)
        totrain=e.rainfall[raindset][:,rainind]
        train = np.nansum(totrain);etrain.append(train)
    edts = np.asarray(edts)
    etrain = np.asarray(etrain)

    # BUILD SEASONAL CYCLES FROM ALL YEARS
    yrs=years
    #if years: yrs=years
    #else: yrs = np.unique(edts[:,0])
    #yrs = np.unique(edts[:,0])
    scycle_train = np.zeros((len(yrs),12))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        for imn in xrange(len(season)):
            mn=season[imn]
            if imn < 5:
                ix = np.where((edts[:,0]==yr) & (edts[:,1]==mn))[0]
            else:
                ix = np.where((edts[:,0]==yr+1) & (edts[:,1]==mn))[0]
            scycle_train[iyr,imn] = np.nansum(etrain[ix])

    return scycle_train, yrs

def scycle_rainsum_raw(rain,dtime,raindset='trmm',season=[8,9,10,11,12,1,2,3,4,5,6,7],years=False):
    '''Calculate seasonal cycle of rainfall from raw dataset (i.e. not only under CBs)
    Sum of all rainfall
    Designed for TRMM rainfall
    rain data in shape (time,lat,lon)
    dtime data in shape (4,ntimesteps)
    RJ July 2016 based on scycle_rainfall
    Can specify years or allow autodiscovery'''

    # SUM OVER LAT AND LONG FOR EVERY TSTEP
    ntime=len(dtime[:,0])
    rfldsum=np.zeros(ntime,dtype=np.float32)
    for t in xrange(ntime):
        rfldsum[t]=np.nansum(rain[t,:,:])

    # BUILD SEASONAL CYCLES FROM ALL YEARS
    yrs=years
    #if years: yrs=years
    #else: yrs = np.unique(dtime[:,0])
    scycle_train = np.zeros((len(yrs),12))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        for imn in xrange(len(season)):
            mn=season[imn]
            if imn < 5:
                ix = np.where((dtime[:,0]==yr) & (dtime[:,1]==mn))[0]
            else:
                ix = np.where((dtime[:,0]==yr+1) & (dtime[:,1]==mn))[0]
            scycle_train[iyr,imn] = np.nansum(rfldsum[ix])

    return scycle_train, yrs


def seasmean(rain, dtime, seas='ann', years=False):
    '''Calculate seasonal mean rainfall from raw dataset
    rain data in shape (time,lat,lon)
    dtime data in shape (4,ntimesteps)
    Can specify years or allow autodiscovery'''

    # Get season
    if seas == 'ann': season=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    if seas == 'djf': season=[12, 1, 2]
    if seas == 'jja': season=[6, 7, 8]

    # MEAN OVER LAT AND LONG FOR EVERY TSTEP
    ntime = len(dtime[:, 0])
    rfldmean = np.zeros(ntime, dtype=np.float32)
    for t in xrange(ntime):
        rfldmean[t] = np.nanmean(rain[t, :, :])

    # SELECT TSTEPS YOU WANT
    if years:
        yrs = years
    else:
        yrs = np.unique(dtime[:, 0])
    seassel = np.zeros((len(season) * len(yrs)),dtype=np.float32)
    z = 0
    for iyr in xrange(len(yrs)):
        yr = yrs[iyr]
        for imn in xrange(len(season)):
            mn = season[imn]
            ix = np.where((dtime[:, 0] == yr) & (dtime[:, 1] == mn))[0]
            seassel[z] = np.nanmean(rfldmean[ix])
            z += 1

    seasmean = np.nanmean(seassel)

    return seasmean

def plotallseasons(scycle,yrs,type='pcolor',anomaly=False,srainfall=False,descr='blank stare'):
    '''Type can be line or pcolor, very different results
    pcolor assumes given a years by months grid'''
    colims=(0,10)
    cm=plt.cm.RdBu
    if anomaly:
        ascycle=scycle-np.tile(scycle.mean(0)[np.newaxis,:],(len(yrs),1))
        colims=(-4,4.5)
        cm=plt.cm.RdBu
        #cm=plt.cm.bwr_r  # Can only use this in matplotlib version 1.0.0 and higher
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']
    if type=='line':
        zz='nothing'
        plt.figure()
        for iyr in xrange(len(yrs)):
            plt.plot(np.arange(1,13),scycle[iyr,:])
            plt.xticks(np.arange(1,13),monthstr)
            plt.ylim(0,6)
            plt.title('Cloudband Season:'+str(yrs[iyr])[2:]+"/"+str(yrs[iyr]+1)[2:])
            plt.grid('on')
            plt.show()
            raw_input()
            #plt.clf()
    if type=='pcolor':
        fig=plt.figure(figsize=[10,6])
        y=np.arange(1,13)[2:9]
        x=yrs
        yy=np.arange(1,14)[2:9]
        xx=np.append(yrs,2011)
        yy,xx = np.meshgrid(yy,xx)
        zz=ascycle[:,2:9]
        ax1 = fig.add_axes([0.05,0.1,0.85,0.5])
        ax2 = fig.add_axes([0.05,0.62,0.85,0.3])
        plt.title('Interannual: '+descr.upper(), fontweight='demibold')
        axcol = fig.add_axes([0.92,0.1,0.03,0.5])
        plt.axes(ax1)
        #plt.contourf(x,y,zz.T,cmap=cm);plt.colorbar(cax=axcol);plt.clim(colims)
        plt.pcolor(xx,yy,zz,cmap=cm);plt.clim(colims)
        plt.yticks(np.arange(1,13)[2:9]+0.5,monthstr[2:9])
        plt.xticks(yrs,yrs,rotation=70)
        plt.grid()
        plt.axis('tight')
        bounds=np.arange(-4.5,5)
        vals=np.arange(-4,5)
        tks=list(vals)
        plt.colorbar(cax=axcol,boundaries=bounds,values=vals,ticks=vals)
        plt.axes(axcol)
        plt.ylabel('anomalous # cloudbands',fontweight='demibold')

        plt.axes(ax2)
        cnts=['0.8','r']
        if srainfall:
            cnt=0
            for tup in srainfall:
                axtra = plt.twinx()
                ryrs, srain = tup
                axtra.bar(ryrs+0.5,srain,width=0.5,align='edge',color=cnts[cnt],alpha=0.4)
                axtra.set_ylim(0,500)
                plt.axes(axtra)
                #plt.yticks([0,250,500,750,1000])
                plt.yticks([0,50,100,150,200,250])
                axtra.set_ylabel('Total Rainfall mm')
                plt.axes(ax2)
                cnt+=1
        n=scycle[:,2:9].sum(1)
        plt.plot(yrs+0.5,n,'b-',lw=1)
        plt.xticks(yrs,[])
        plt.xlim(yrs[0],yrs[-1]+1)
        plt.ylim(5,25)
        plt.grid()

        fname=figdir+'/2panelSeason-'+descr.strip()+'.png'
        #plt.savefig(fname,dpi=200)
        #plt.axis('tight')
        #plt.show()
    return zz

def plotallseasonsRain(scycle,yrs,type='pcolor',anomaly=False,srainfall=False,\
    descr='blank stare'):
    '''Type can be line or pcolor, very different results
    pcolor assumes given a years by months grid'''
    # alternative line option added by RJ 'line_rj'
    # I think I added this because I couldn't get 'line' to work!
    colims=(0,10)
    cm=plt.cm.RdBu
    if anomaly:
        ascycle=scycle-np.tile(scycle.mean(0)[np.newaxis,:],(len(yrs),1))
        #colims=(-200,200)
        cm=plt.cm.RdBu
        #cm=plt.cm.bwr_r
    fig=plt.figure(figsize=[10,6])
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar',\
              'Apr','May','Jun','Jul']
    if type=='line':
        for iyr in xrange(len(yrs)):
            plt.plot(np.arange(1,13),scycle[iyr,:])
            plt.xticks(np.arange(1,13),monthstr)
            #plt.ylim(0,10)
            plt.title('Cloudband Season:'+str(yrs[iyr])[2:]+"/"+str(yrs[iyr]+1)[2:])
            plt.grid('on')
            plt.show()
            raw_input()
            plt.clf()
    if type=='line_rj':
        ymonsum = np.zeros(12, dtype=np.float32)
        for imn in xrange(len(monthstr)):
            ymonsum[imn] = np.nansum(scycle[:, imn])
        plt.plot(np.arange(0, 12), ymonsum)
        plt.xticks(np.arange(0, 12), monthstr)
        plt.xlim(0, 11)
        # plt.ylim(0,10)
        plt.grid('on')
        fname = 'SeasonalCycle_totalTTCBrain_' + descr.strip() + '.png'
        plt.savefig(fname, dpi=200)
    if type=='pcolor':
        y=np.arange(1,13)[2:9]
        x=yrs
        yy=np.arange(1,14)[2:9]
        xx=np.append(yrs,2011)
        yy,xx = np.meshgrid(yy,xx)
        zz=ascycle[:,2:9]
        ax1 = fig.add_axes([0.05,0.1,0.85,0.5])
        ax2 = fig.add_axes([0.05,0.62,0.85,0.3])
        plt.title('Interannual: '+descr.upper(), fontweight='demibold')
        axcol = fig.add_axes([0.92,0.1,0.03,0.5])
        plt.axes(ax1)
        #plt.contourf(x,y,zz.T,cmap=cm);plt.colorbar(cax=axcol);plt.clim(colims)
        plt.pcolor(xx,yy,zz,cmap=cm);plt.clim(colims)
        plt.yticks(np.arange(1,13)[2:9]+0.5,monthstr[2:9])
        plt.xticks(yrs,yrs,rotation=70)
        plt.grid()
        plt.axis('tight')
        bounds=np.arange(-4.5,5)
        vals=np.arange(-4,5)
        tks=list(vals)
        #plt.colorbar(cax=axcol,boundaries=bounds,values=vals,ticks=vals)
        plt.colorbar(cax=axcol)
        plt.axes(axcol)
        plt.ylabel('anomalous # cloudbands',fontweight='demibold')

        plt.axes(ax2)
        cnts=['0.8','r']
        if srainfall:
            cnt=0
            for tup in srainfall:
                axtra = plt.twinx()
                ryrs, srain = tup
                axtra.bar(ryrs+0.5,srain,width=0.5,align='edge',color=cnts[cnt],alpha=0.4)
                axtra.set_ylim(0,500)
                plt.axes(axtra)
                #plt.yticks([0,250,500,750,1000])
                plt.yticks([0,50,100,150,200,250])
                axtra.set_ylabel('Total Rainfall mm')
                plt.axes(ax2)
                cnt+=1
        n=scycle[:,2:9].sum(1)
        plt.plot(yrs+0.5,n,'b-',lw=1)
        plt.xticks(yrs,[])
        plt.xlim(yrs[0],yrs[-1]+1)
        plt.ylim(5,25)
        plt.grid()

        fname=figdir+'/2panelSeason-'+descr.strip()+'.png'
        #plt.savefig(fname,dpi=200)
        #plt.axis('tight')
        #plt.show()
    return zz



def plotseasonbox(scycle,descr,ax=False,savefig=False,ylims=False):
    if isinstance(ax,bool):plt.figure()
    else: plt.axes(ax)
    lc='b'
    bp = dict(linewidth=2,color=lc)
    fp = dict(markersize=12,markeredgewidth=1.5,markeredgecolor=lc,alpha=1.)
    mp = dict(linewidth=2,color='r')
    #wp = dict(linewidth=1.5,color=lc,alpha=alpha)
    wp = dict(linewidth=1.5,color=lc)
    cp = dict(linewidth=1.5,color=lc)
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb',\
              'Mar','Apr','May','Jun','Jul']
    plt.boxplot(scycle, notch=0, sym='+', vert=1, whis=1.5,\
                showfliers=True,showcaps=True,boxprops=bp,flierprops=fp,\
                medianprops=mp,whiskerprops=wp,capprops=cp)
    #plt.violinplot(scycle,
    #               showmeans=False,
    #               showmedians=True)
    #plt.plot(np.arange(1,13),np.median(scycle,axis=0),'k-',lw=1.5)
    plt.xticks(np.arange(1,13),monthstr,fontsize=14.0)
    plt.yticks(np.arange(1,15),fontsize=14.0)
    if ylims:
        plt.ylim(ylims[0],ylims[1])
    else:
        plt.ylim(0,10.)
    plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
    plt.title(descr, fontweight='demibold')
    plt.grid()
    fname=figdir+'scycle-'+descr+'.png'
    if savefig: plt.savefig(fname,dpi=150)
    plt.show()

def plotseasonbox_background(scycle,ax=False,savefig=False,ylims=False):
    if isinstance(ax,bool):plt.figure()
    else: plt.axes(ax)
    bp = dict(linewidth=2,color='None',facecolor='.8')
    mp = dict(linewidth=2,color='.5')
    wp = dict(color='None')
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb',\
              'Mar','Apr','May','Jun','Jul']
    plt.boxplot(scycle, notch=0, widths=.4, vert=1,patch_artist=True,\
                showfliers=False,showcaps=False,boxprops=bp,\
                medianprops=mp,whiskerprops=wp)
    #plt.violinplot(scycle,
    #               showmeans=False,
    #               showmedians=True)
    #plt.plot(np.arange(1,13),np.median(scycle,axis=0),'k-',lw=1.5)
    plt.xticks(np.arange(1,13),monthstr,fontsize=14.0)
    plt.yticks(np.arange(1,15),fontsize=14.0)
    if ylims:
        plt.ylim(ylims[0],ylims[1])
    else:
        plt.ylim(0,10.)

def plotseasonbox_rj(scycle,descr,picext,savefig=False):
    # Alternative plotseasonbox made by RJ
    # descr is title of the plot
    # pixext is the beginning of file name (can include dir)
    plt.figure()
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']
    plt.boxplot(scycle, notch=0, sym='+', vert=1, whis=1.5) # produces boxplot
    plt.plot(np.arange(1,13),scycle.mean(0),'k-',lw=1) # produces mean line
    plt.xticks(np.arange(1,13),monthstr,fontsize=13.0) # month labels
    plt.yticks(np.arange(1,14),fontsize=13.0)
    plt.ylim(0,8.5)
    plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
    plt.title(descr.upper(), fontweight='demibold')
    #plt.grid()
    fname=picext+'_scycle.png'
    if savefig: plt.savefig(fname,dpi=150)

def plotseasonbox_rain(scycle,descr,ax=False):
    if isinstance(ax,bool):plt.figure()
    else: plt.axes(ax)
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']
    plt.boxplot(scycle, notch=0, sym='+', vert=1, whis=1.5)
    plt.plot(np.arange(1,13),st.nanmean(scycle,0),'k-',lw=1)
    plt.xticks(np.arange(1,13),monthstr,fontsize=13.0)
    #plt.yticks(np.arange(1,10),fontsize=13.0)
    #plt.ylim(0,9.5)
    #plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
    plt.title('Seasonal Cycle: '+descr.upper(), fontweight='demibold')
    fname=figdir+'scycle-'+descr+'.png'
    #plt.savefig(fname,dpi=200)
    plt.show()

def plotseasonbox_wrain(ntttcycle,raincycle,descr,dset,pcent=True,savefig=False):
    '''Plots seasonal cycle of number of TTTs and
     rainfall contributed by TTTs'''
    # RJ August 2016

    fig, ax1 = plt.subplots()
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']

    # Plot scycle of number of TTTs
    ax1.boxplot(ntttcycle, notch=0, sym='+', vert=1, whis=1.5) # produces boxplot
    ax1.plot(np.arange(1,13),ntttcycle.mean(0),'k-',lw=1,zorder=3) # produces mean line
    plt.xticks(np.arange(1,13),monthstr,fontsize=13.0) # month labels
    plt.yticks(np.arange(1,14),fontsize=13.0)
    plt.ylim(0,12.5)
    ax1.set_ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
    plt.title(descr.upper(), fontweight='demibold')
    #plt.grid()

    # Plot scycle of rainfall
    ax2 = ax1.twinx()
    ax2.fill_between(np.arange(1,13),0,raincycle,edgecolor='none',facecolor='lightskyblue',zorder=1)
    if pcent:
        ax2.set_ylabel('%',color='black')
        ax2.set_ylim(0,100)
    else:
        ax2.set_ylabel('mm',color='black')
    ax1.set_zorder(2)
    ax2.set_zorder(1)
    ax1.set_axis_bgcolor('none')

    fname='Scycle_nttt_rain-'+descr+'_'+dset+'.png'
    if savefig: plt.savefig(fname,dpi=150)

def plotseasonbox_wrain_2(ntttcycle,raincycle,descr,dset,pcent=True,savefig=False):
    '''Plots seasonal cycle of number of TTTs and
         rainfall contributed by TTTs'''
    # RJ August 2016

    fig, ax1 = plt.subplots()
    monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']

    # Plot scycle of number of TTTs
    boxprops = dict(linestyle='-',linewidth=1.5,color='k')
    whiskerprops = dict(linestyle='--',linewidth=1.5,color='k')
    medianprops = dict(linestyle='-',linewidth=1.5,color='k')
    capprops = dict(linestyle='-',linewidth=1.5,color='k')
    ax1.boxplot(ntttcycle, notch=0, sym='+', vert=1, whis=1.5, boxprops=boxprops, medianprops=medianprops, whiskerprops=whiskerprops, capprops=capprops) # produces boxplot
    #ax1.plot(np.arange(1,13),ntttcycle.mean(0),'k-',lw=1,zorder=3) # produces mean line
    plt.xticks(np.arange(1,13),monthstr,fontsize=13.0) # month labels
    plt.yticks(np.arange(1,14),fontsize=13.0)
    plt.ylim(0,8.5)
    ax1.set_ylabel('No. of TTT events', fontsize=13.0, weight='demibold')
    #plt.title(descr.upper(), fontweight='demibold')
    #plt.grid()

    # Plot scycle of rainfall
    ax2 = ax1.twinx()
    ax2.fill_between(np.arange(1,13),0,raincycle,edgecolor='none',facecolor='mediumaquamarine',zorder=1)
    if pcent:
        ax2.set_ylabel('%',color='black')
        ax2.set_ylim(0,30)
    else:
        ax2.set_ylabel('mm',color='black')
    ax1.set_zorder(2)
    ax2.set_zorder(1)
    ax1.set_axis_bgcolor('none')

    fname='Scycle_nttt_rain-'+descr+'_'+dset+'.png'
    if savefig: plt.savefig(fname,dpi=150)

def spatiofreq(s,eventkeys,descr,plottrk=False,plothex=False,res=4.0,sub='SA'):
    '''Get grid-cell frequencies for cloudband tracks'''
    plt.close('all')
    mbskeys = s.mbskeys
    refkey = mbskeys[0]
    basekey = mbskeys[1]
    try:
        dset, varstr, levsel, deriv, expid = basekey.split('-')
    except:
        dset, varstr, levsel, deriv = basekey.split('-')
    vkey='%s-%s-%s-%s' %(dset, varstr, levsel, deriv)
    x1,x2,y1,y2=blb.blobfilters[sub+'cloudband'][vkey]['ROI']
    nx, ny = np.abs(x1-x2)/res, np.abs(y1-y2)/res
    grdsz = (np.int32(nx),np.int32(ny))

    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    x,y=np.array(()),np.array(())

    for k in eventkeys:
        e = s.events[k]
        if plottrk: plt.plot(e.trkcX,e.trkcY,'k-x',lw=1,alpha=0.1)
        #x,y = np.append(x,e.trkcX.mean()), np.append(y,e.trkcY.mean())
        x,y = np.append(x,e.trkcX), np.append(y,e.trkcY)

    if plothex:
        cm=plt.cm.bone_r
        plt.hexbin(x,y,gridsize=grdsz,mincnt=1,cmap=cm)
        plt.colorbar()

    plt.xlim(x1,x2)
    plt.ylim(vlt)
    plt.clim(10,80)
    plt.title('Event Spatial Count: '+descr.upper(), fontweight='demibold')
    fname=figdir+'hexbin-'+descr+'.png'
    #plt.savefig(fname,dpi=200)
    #plt.show()

def spatiofreq2(m,s,lat,lon,yrs,eventkeys,meanmask=False,figno=1,\
                clim=(4,36,4),month=False,flagonly=False,fontdict=False):
    '''Get grid-cell frequencies for no. of times a grid-cell falls within a 
       contour describing a feature from metblobs.
    USAGE: If wish to create Basemap within function, m will be "create"
           if wish to have only for particular month, month=yourchoice
           if wish to only count for flagged days, flagonly=True'''
    #plt.close('all')
    if not fontdict: fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbskeys = s.mbskeys
    refkey = s.events.values()[0].refkey
    basekey = refkey
    try:
        dset, varstr, levsel, deriv, expid = basekey.split('-')
        descr = dset+'-'+varstr
    except:
        dset, varstr, levsel, deriv = basekey.split('-')
        descr = dset+'-'+varstr
    #vkey='%s-%s-%s-%s' %(dset, varstr, levsel, deriv)
    #x1,x2,y1,y2=blb.blobfilters[sub+'cloudband'][vkey]['ROI']
    #nx, ny = np.abs(x1-x2)/res, np.abs(y1-y2)/res
    #grdsz = (np.int32(nx),np.int32(ny))

    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    x,y=np.array(()),np.array(())

    #lon=np.arange(lon[0],lon[-1]+1,0.5)
    #lat=np.arange(lat[0],lat[-1]-1,-0.5)
    allmask = np.zeros((lat.shape[0],lon.shape[0]),dtype=np.float32)

    #Determine which variable we dealing with
    countkey=refkey
    if flagonly: countkey = s.flagkey

    for k in eventkeys:
        e = s.events[k]
        if month:
            mn = month
            mst = e.trkdtimes[0,1]
            if mst != mn: continue
        if flagonly:
            itrk=e.ixflags
        else:
            trkarr = np.int32(e.trkarrs[countkey])
            if trkarr.ndim==2:
                ixt = np.where(trkarr[:,1]>0)[0]
                uni,iu=np.unique(trkarr[ixt,0],return_index=True)
                itrk = trkarr[ixt,1]
                #print len(itrk),":",itrk
            elif trkarr.ndim==3:
                itrk = np.ndarray((0,),dtype=np.int32)
                for d in xrange(trkarr.shape[2]):
                    ixt = np.where(trkarr[:,1,d]>0)[0]
                    uni,iu=np.unique(trkarr[ixt,0],return_index=True)
                    ixt=ixt[iu]
                    itrk = np.append(itrk,trkarr[ixt,1,d].squeeze())
        # Get masks for each contour feature of a track
        for ixtrk in itrk:
            mask = my.poly2mask(lon,lat,e.blobs[countkey]['ch'][ixtrk])
            allmask=allmask+np.float32(mask)

    #cm=plt.cm.PuBu
    #cm=plt.cm.gist_earth_r
    #cm=plt.cm.YlGnBu
    #cm=plt.cm.binary
    #cm=plt.cm.OrRd
    if isinstance(meanmask,np.ndarray):cm=plt.cm.RdBu;#cm=plt.cm.bwr
    else:cm=plt.cm.gist_gray_r
    #cm=plt.cm.jet
    #cm=plt.cm.gray_r
    cm.set_under(color='w')
    if m=='create':
        m, f = blb.SAfrBasemap(lat[3:-6],lon[3:-3],drawstuff=True,\
                               prj='cyl',fno=figno,rsltn='l')
    #df=m.transform_scalar(allmask[::-1,:],lon,lat[::-1],len(lon),len(lat))
    #m.imshow(df,cm,interpolation='nearest')
    #m.pcolor(lon,lat,allmask,cmap=cm)
    #plt.clim(300,1200)
    lon,lat = np.meshgrid(lon,lat)
    #m.contourf(lon,lat,(allmask/len(yrs)),cmap=cm)
    #
    ### Next is dirty trick to drop out permanent mid-latitude cloudiness
    ### allowing colormap to enhance cloud bands better
    allmask[-9:,:]=0 
    std_mask=allmask/len(yrs)
    if isinstance(m,bool):
        return std_mask
    # IF M WAS NOT A BOOLEAN (False) THIS WILL CONTINUE TO PLOTTING
    if isinstance(meanmask,np.ndarray):
        std_mask=std_mask-meanmask
        std_mask=np.where(np.abs(std_mask)<.5,np.nan,std_mask)
        lnmin,lnmx,latmn,latmx =\
                        blb.filters.blobfilters['SAcloudband'][countkey]['ROI']
        latmask = (lat<latmn) & (lat>latmx) # this is for S. Hemisphere
        meanmasked = np.ma.MaskedArray(meanmask,mask=~latmask)
        m.contour(lon,lat,meanmasked,[2,4],colors='k')
        #m.contourf(lon,lat,meanmasked,[2,4,14],hatches=['.','..'],\
        #          colors='none',linestyles='-',linewidths='1',alpha=.1)
    ## NEED TO DO THIS SINCE PCOLOR IS NOT SHADING VALUES OUTSIDE OF THE CLIMS
    cstd_mask=np.where(std_mask>clim[1],clim[1],std_mask)
    cstd_mask=np.where(cstd_mask<clim[0],clim[0],cstd_mask)
    # Plot pcolor
    pcolmap=m.pcolormesh(lon,lat,cstd_mask,cmap=cm,zorder=1)
    img=plt.gci()

    #Plotting centroids
    for k in eventkeys:
        e = s.events[k]
        if month:
            mn = month
            mst = e.trkdtimes[0,1]
            if mst != mn: continue
        m.plot(e.trkcX[0],e.trkcY[0],color='w',marker='o',markersize=4)
        if 'COL' in e.mbskeys:
            if len(e.assoctrks['COL'])>0:
                trar=e.trkarrs['COL']
                for ntrk in xrange(trar.shape[2]):
                    ixx = np.where(trar[:,1,ntrk]>0)[0]
                    xx, yy = trar[ixx,2,ntrk],trar[ixx,3,ntrk]
                    off=np.random.rand()*.5
                    m.scatter(xx[-1]+off,yy[-1]+off,50,color='none',\
                              edgecolor='k',marker='o',linewidth='2')
                    plt.sci(img)
        #m.plot(e.trkcX,e.trkcY,'0.5')

    plt.clim(clim[0],clim[1])
    bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    #vals=np.arange(0,35,2)
    if not month:
        f,ax=plt.gcf(),plt.gca()
        axcol=f.add_axes([0.93,0.2,0.02,0.6])
        plt.colorbar(mappable=img,cax=axcol,boundaries=bounds)
        my.ytickfonts()
        if isinstance(meanmask,np.ndarray):
            plt.ylabel('anomaly grid-point count / year',fontdict=fd)
        else:
            plt.ylabel('grid-point count / year',fontdict=fd)
        plt.axes(ax)
        plt.title('Cloudband Annual Grid-Point Count Climatology: '\
                  +descr.upper(),fontsize='14',fontdict=fd)
        fname=figdir+'/FootprintFreqencygray-'+descr+'.png'
        if flagonly:
            fname=figdir+'/FootprintFreqencygray-'+descr+'_flagonly.png'
        plt.savefig(fname,dpi=150)
    elif month:
        f,ax=plt.gcf(),plt.gca()
        axcol=f.add_axes([0.93,0.2,0.02,0.6])
        plt.colorbar(cax=axcol,boundaries=bounds)
        my.ytickfonts()
        if isinstance(meanmask,np.ndarray):
            plt.ylabel('anomaly grid-point count / year',fontdict=fd)
        else:
            plt.ylabel('grid-point count / year',fontdict=fd)
        plt.axes(ax)
        plt.title(mndict[month], fontweight='demibold')

    return std_mask

def spatiofreq3(m,s,lat,lon,yrs,eventkeys,meanmask=False,figno=1,\
                clim=(4,36,4),month=False,flagonly=False,fontdict=False):
    '''Get grid-cell frequencies for no. of times a grid-cell falls within a
       contour describing a feature from metblobs.
       spatiofreq3 similar to spatiofreq2 but edited by RJ
    USAGE: If wish to create Basemap within function, m will be "create"
           if wish to have only for particular month, month=yourchoice
           if wish to only count for flagged days, flagonly=True'''
    if not fontdict: fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbskeys = s.mbskeys
    refkey = s.events.values()[0].refkey
    basekey = refkey
    try:
        dset, varstr, levsel, deriv, expid = basekey.split('-')
        descr = dset+'-'+varstr
    except:
        dset, varstr, levsel, deriv = basekey.split('-')
        descr = dset+'-'+varstr
    #vkey='%s-%s-%s-%s' %(dset, varstr, levsel, deriv)
    #x1,x2,y1,y2=blb.blobfilters[sub+'cloudband'][vkey]['ROI']
    #nx, ny = np.abs(x1-x2)/res, np.abs(y1-y2)/res
    #grdsz = (np.int32(nx),np.int32(ny))

    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    x,y=np.array(()),np.array(())

    #lon=np.arange(lon[0],lon[-1]+1,0.5)
    #lat=np.arange(lat[0],lat[-1]-1,-0.5)
    allmask = np.zeros((lat.shape[0],lon.shape[0]),dtype=np.float32)

    #Determine which variable we dealing with
    countkey=refkey
    if flagonly: countkey = s.flagkey

    for k in eventkeys:
        e = s.events[k]
        if month:
            mn = month
            mst = e.trkdtimes[0,1]
            if mst != mn: continue
        if flagonly:
            itrk=e.ixflags
        else:
            trkarr = np.int32(e.trkarrs[countkey])
            if trkarr.ndim==2:
                ixt = np.where(trkarr[:,1]>0)[0]
                uni,iu=np.unique(trkarr[ixt,0],return_index=True)
                itrk = trkarr[ixt,1]
                #print len(itrk),":",itrk
            elif trkarr.ndim==3:
                itrk = np.ndarray((0,),dtype=np.int32)
                for d in xrange(trkarr.shape[2]):
                    ixt = np.where(trkarr[:,1,d]>0)[0]
                    uni,iu=np.unique(trkarr[ixt,0],return_index=True)
                    ixt=ixt[iu]
                    itrk = np.append(itrk,trkarr[ixt,1,d].squeeze())
        # Get masks for each contour feature of a track
        for ixtrk in itrk:
            mask = my.poly2mask(lon,lat,e.blobs[countkey]['ch'][ixtrk])
            allmask=allmask+np.float32(mask)

    #cm=plt.cm.PuBu
    #cm=plt.cm.gist_earth_r
    #cm=plt.cm.YlGnBu
    #cm=plt.cm.binary
    #cm=plt.cm.OrRd
    if isinstance(meanmask,np.ndarray):cm=plt.cm.RdBu;#cm=plt.cm.bwr
    else:cm=plt.cm.gist_gray_r
    #cm=plt.cm.jet
    #cm=plt.cm.gray_r
    cm.set_under(color='w')
    if m=='create':
        m, f = blb.SAfrBasemap(lat[3:-6],lon[3:-3],drawstuff=True,\
                               prj='cyl',fno=figno,rsltn='l')
    #df=m.transform_scalar(allmask[::-1,:],lon,lat[::-1],len(lon),len(lat))
    #m.imshow(df,cm,interpolation='nearest')
    #m.pcolor(lon,lat,allmask,cmap=cm)
    #plt.clim(300,1200)
    lon,lat = np.meshgrid(lon,lat)
    #m.contourf(lon,lat,(allmask/len(yrs)),cmap=cm)
    #
    ### Next is dirty trick to drop out permanent mid-latitude cloudiness
    ### allowing colormap to enhance cloud bands better
    allmask[-9:,:]=0 # the -9: sets all of the values from 9 below the bottom to 0
                        # for the SA domain and NOAA res this is to 37.5
                        # I dont think it does anything because these parts are excluded anyway!
    std_mask=allmask/len(yrs)
    if isinstance(m,bool):
        return std_mask
    # IF M WAS NOT A BOOLEAN (False) THIS WILL CONTINUE TO PLOTTING
    if isinstance(meanmask,np.ndarray):
        std_mask=std_mask-meanmask
        std_mask=np.where(np.abs(std_mask)<.5,np.nan,std_mask)
        lnmin,lnmx,latmn,latmx =\
                        blb.filters.blobfilters['SAcloudband'][countkey]['ROI']
        latmask = (lat<latmn) & (lat>latmx) # this is for S. Hemisphere
        meanmasked = np.ma.MaskedArray(meanmask,mask=~latmask)
        m.contour(lon,lat,meanmasked,[2,4],colors='k')
        #m.contourf(lon,lat,meanmasked,[2,4,14],hatches=['.','..'],\
        #          colors='none',linestyles='-',linewidths='1',alpha=.1)
    ## NEED TO DO THIS SINCE PCOLOR IS NOT SHADING VALUES OUTSIDE OF THE CLIMS
    cstd_mask=np.where(std_mask>clim[1],clim[1],std_mask)
    cstd_mask=np.where(cstd_mask<clim[0],clim[0],cstd_mask)
    # Plot pcolor
    pcolmap=m.pcolormesh(lon,lat,cstd_mask,cmap=cm,zorder=1)
    img=plt.gci()

    #Plotting centroids
    # for k in eventkeys:
    #     e = s.events[k]
    #     if month:
    #         mn = month
    #         mst = e.trkdtimes[0,1]
    #         if mst != mn: continue
    #     m.plot(e.trkcX[0],e.trkcY[0],color='w',marker='o',markersize=4)
    #     if 'COL' in e.mbskeys:
    #         if len(e.assoctrks['COL'])>0:
    #             trar=e.trkarrs['COL']
    #             for ntrk in xrange(trar.shape[2]):
    #                 ixx = np.where(trar[:,1,ntrk]>0)[0]
    #                 xx, yy = trar[ixx,2,ntrk],trar[ixx,3,ntrk]
    #                 off=np.random.rand()*.5
    #                 m.scatter(xx[-1]+off,yy[-1]+off,50,color='none',\
    #                           edgecolor='k',marker='o',linewidth='2')
    #                 plt.sci(img)
    #     #m.plot(e.trkcX,e.trkcY,'0.5')

    plt.clim(clim[0],clim[1])
    bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    #vals=np.arange(0,35,2)
    if not month:
        f,ax=plt.gcf(),plt.gca()
        axcol=f.add_axes([0.93,0.2,0.02,0.6])
        plt.colorbar(mappable=img,cax=axcol,boundaries=bounds)
        my.ytickfonts()
        if isinstance(meanmask,np.ndarray):
            plt.ylabel('anomaly grid-point count / year',fontdict=fd)
        else:
            plt.ylabel('grid-point count / year',fontdict=fd)
        plt.axes(ax)
        plt.title('Cloudband Annual Grid-Point Count Climatology: '\
                  +descr.upper(),fontsize='14',fontdict=fd)
        fname=figdir+'/FootprintFreqencygray-'+descr+'.png'
        if flagonly:
            fname=figdir+'/FootprintFreqencygray-'+descr+'_flagonly.png'
        plt.savefig(fname,dpi=150)
    elif month:
        f,ax=plt.gcf(),plt.gca()
        axcol=f.add_axes([0.93,0.2,0.02,0.6])
        plt.colorbar(cax=axcol,boundaries=bounds)
        my.ytickfonts()
        if isinstance(meanmask,np.ndarray):
            plt.ylabel('anomaly grid-point count / year',fontdict=fd)
        else:
            plt.ylabel('grid-point count / year',fontdict=fd)
        plt.axes(ax)
        plt.title(mndict[month], fontweight='demibold')

    return std_mask

def spatiofreq4(m,s,modname,lat,lon,yrs,eventkeys,per='year',meanmask=False,\
                clim=(4,36,4),month=False,savefig=False,\
                flagonly=False,col='col',cens=False,frm_event='all'):
    '''Get grid-cell frequencies for no. of times a grid-cell falls within a
       contour describing a feature from metblobs.
       spatiofreq4 similar to spatiofreq2 but edited by RJ for CMIP5 multiplotting
       per is used to determine if it's plotting per year or per CBs in this model
       cens is for plotting centroids - if None, no cens, if 'All' all, if array, that array
    USAGE: if wish to have only for particular month, month=yourchoice
           if wish to only count for flagged days, flagonly=True'''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    allmask = np.zeros((lat.shape[0],lon.shape[0]),dtype=np.float32)

    #Determine which variable we dealing with
    countkey=s.events.values()[0].refkey
    if flagonly: countkey = s.flagkey

    cnt=0
    for k in eventkeys:
        e = s.events[k]
        if month:
            mn = month
            mst = e.trkdtimes[0,1]
            if mst != mn: continue
        if flagonly:
            itrk=e.ixflags
        else:
            trkarr = np.int32(e.trkarrs[countkey])
            if trkarr.ndim==2:
                ixt = np.where(trkarr[:,1]>0)[0]
                itrk = trkarr[ixt,1]
            elif trkarr.ndim==3:
                itrk = np.ndarray((0,),dtype=np.int32)
                for d in xrange(trkarr.shape[2]):
                    ixt = np.where(trkarr[:,1,d]>0)[0]
                    uni,iu=np.unique(trkarr[ixt,0],return_index=True)
                    ixt=ixt[iu]
                    itrk = np.append(itrk,trkarr[ixt,1,d].squeeze())
        # Get masks for each contour feature of a track
        if frm_event=='all':
            for ixtrk in itrk:
                mask = my.poly2mask(lon,lat,e.blobs[countkey]['ch'][ixtrk])
                allmask=allmask+np.float32(mask)
                cnt+=1
        elif frm_event=='first':
            ixtrk=itrk[0]
            mask = my.poly2mask(lon, lat, e.blobs[countkey]['ch'][ixtrk])
            allmask = allmask + np.float32(mask)
            cnt += 1

    nblobs=cnt

    if isinstance(meanmask,np.ndarray):cm=plt.cm.RdBu;
    else:
        if col=='col':
            cm = plt.cm.magma
        elif col=='bw':
            cm=plt.cm.gist_gray_r
    #cm.set_under(color='w')

    if per=='year':
        std_mask=allmask/len(yrs)
    elif per=='cbs':
        std_mask=allmask/nblobs*100
    if isinstance(meanmask,np.ndarray):
        std_mask=std_mask-meanmask
        std_mask=np.where(np.abs(std_mask)<.5,np.nan,std_mask)
        lnmin,lnmx,latmn,latmx =\
                        blb.filters.blobfilters['SAcloudband'][countkey]['ROI']
        latmask = (lat<latmn) & (lat>latmx) # this is for S. Hemisphere
        meanmasked = np.ma.MaskedArray(meanmask,mask=~latmask)
        m.contour(lon,lat,meanmasked,[2,4],colors='k')

    ## NEED TO DO THIS SINCE PCOLOR IS NOT SHADING VALUES OUTSIDE OF THE CLIMS
    cstd_mask=np.where(std_mask>clim[1],clim[1],std_mask)
    cstd_mask=np.where(cstd_mask<clim[0],clim[0],cstd_mask)
    # Plot pcolor
    pcolmap=m.pcolormesh(lon,lat,cstd_mask,cmap=cm,zorder=1)
    img=plt.gci() # gets a reference for the image

    if cens !='None':
        #Plotting centroids
        if cens=='all':
            print 'Plotting all centroids'
            for k in eventkeys:
                e = s.events[k]
                if month:
                    mn = month
                    mst = e.trkdtimes[0,1]
                    if mst != mn: continue
                m.plot(e.trkcX[0],e.trkcY[0],color='k',marker='o',markersize=0.5)
        else:
            print 'Plotting all centroids'
            for k in eventkeys:
                e = s.events[k]
                if month:
                    mn = month
                    mst = e.trkdtimes[0,1]
                    if mst != mn: continue
                m.plot(e.trkcX[0],e.trkcY[0],color='k',marker='o',markersize=0.5)

            print 'Plotting centroids for sample '
            cont_cens=cens[0]
            mada_cens=cens[1]
            m.plot(cont_cens[:,0], cont_cens[:,1], color='fuchsia', marker='o', markersize=0.5,markeredgecolor='fuchsia',linestyle='None')
            m.plot(mada_cens[:,0], mada_cens[:,1], color='blue', marker='o', markersize=0.5,markeredgecolor='blue',linestyle='None')


    plt.clim(clim[0],clim[1]) # sets color limits of current image
    bounds=np.arange(clim[0],clim[1]+clim[2],clim[2])
    if not month:
        if savefig:
            f,ax=plt.gcf(),plt.gca()
            axcol=f.add_axes([0.93,0.2,0.02,0.6])
            plt.colorbar(mappable=img,cax=axcol,boundaries=bounds)
            my.ytickfonts()
            if isinstance(meanmask,np.ndarray):
                plt.ylabel('anomaly grid-point count / year',fontdict=fd)
            else:
                plt.ylabel('grid-point count / year',fontdict=fd)
            plt.axes(ax)
            plt.title('Cloudband Grid-Point Count Climatology: '\
                      +descr.upper(),fontsize='14',fontdict=fd)
            fname='/FootprintFreqencygray-'+descr+'.png'
            if flagonly:
                fname='/FootprintFreqencygray-'+descr+'_flagonly.png'
            plt.savefig(fname,dpi=150)
        else:
            f, ax = plt.gcf(), plt.gca() # get reference and set axes
            axcol = f.add_axes([0.91, 0.15, 0.01, 0.6])
            plt.colorbar(cax=axcol, boundaries=bounds)
            my.ytickfonts(fontweight='normal',fontsize=10)
            if isinstance(meanmask, np.ndarray):
                plt.ylabel('anomaly grid-point count / year', fontsize=10)
            else:
                if per=='year':
                    plt.ylabel('grid-point count / year', fontsize=10)
                elif per=='cbs':
                    plt.ylabel('% of cbs covering gridbox', fontsize=10)
            plt.axes(ax)
            plt.title(modname,fontsize=8)
    elif month:
        f,ax=plt.gcf(),plt.gca()
        axcol=f.add_axes([0.93,0.2,0.02,0.6])
        plt.colorbar(cax=axcol,boundaries=bounds)
        my.ytickfonts()
        if isinstance(meanmask,np.ndarray):
            plt.ylabel('anomaly grid-point count / year',fontsize=10)
        else:
            plt.ylabel('grid-point count / year',fontsize=10)
        plt.axes(ax)
        plt.title(mndict[month], fontweight='demibold')

    return std_mask


def wetlandevents(s, eventkeys, land_datasets=['wrc']):
    '''Returns the keys of events that produce continental rainfall'''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    ### Select events based on the criteria of continental rainfall 
    ### using SA datasets
    wetlandkeys=[]
    for k in eventkeys:
        e = s.events[k]
        for dta in land_datasets:
            rain = e.rainfall[dta]
            if (rain.shape > 0) & np.any( ~np.isnan(rain)):
                #print "Wet event in",dta," on:",e.trkdtimes[0]
                wetlandkeys.append(k)
    wetlandkeys.sort()
    wetland = np.unique(np.asarray(wetlandkeys))
    wetlandkeys = list(wetland)

    return wetlandkeys

def blobcontours(s,eventkeys,chkey='noaa-olr-0-all'):
    '''Simply returns a dictionary of a list of polygons (ndarrays) describing outline of chosen variable
    for each event
    Takes this directly out of blobs dictionary eg blobs[blobkey]['ch'][ix]
    to return something a little more manageable'''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])
    allpolys = {}
    for k in eventkeys:
        e=s.events[k]
        #trkarr = np.int32(e.trkarrs[chkey][])
        polys = []
        for ix in itrk:
            polys.append(s.blobs[chkey]['ch'][ix])
        allpolys[k] = polys

    return allpolys

def raineventmask(allpolys,s,raindata):
    '''Returns a mask for rainfall for each day of each event
    format: dictionary for each entry
    Currently only uses WRC stations'''
    keylist = allpolys.keys()
    keylist.sort()

    rain,date,xypts = raindata
    rainmasks = {}
    masklist=[]
    maskdtimes=[]
    maskekeys=[]
    for k in keylist:
        polylist=[]
        dtimes=s.events[k].trkdtimes
        cnt=0
        for poly in allpolys[k]:
            chmask = points_inside_poly(xypts,poly)
            polylist.append(chmask)
            masklist.append(chmask)
            maskdtimes.append(dtimes[cnt])
            maskekeys.append(k)
            cnt+=1
        rainmasks[k] = polylist

    return rainmasks, masklist, maskdtimes, maskekeys

# def griddedrainmasks(s,eventkeys,raindata,refkey='noaa-olr-0-0'):
#     '''Returns a masked array of rainfall from TTT events (tstep, lon, lat)'''
#     # RJ 2016
#     # THIS IS BROKEN _ MASKED ARRAY FAILS - USED ALTERNATIVE APPROACH IN AP
#     #input raindata should be in the form raindata=(rain,dtime,(lon,lat))
#     #if want to choose particular years or mons need to do this with eventkey input
#     #e.g. using stats.specificseason
#     if not eventkeys:
#         eventkeys=[]
#         for ed in s.uniques:
#             eventkeys.append(ed[0])
#
#     rain,date,xypts = raindata
#     lon, lat = xypts
#     nlon=len(lon)
#     nlat=len(lat)
#     routput=[]
#     for k in eventkeys:
#         e=s.events[k]
#         rtmp=np.ma.zeros(((len(e.trkdtimes)),nlat,nlon),dtype=np.float32)
#         for t in xrange(len(e.trkdtimes)):
#             ix = my.ixdtimes(date,[e.trkdtimes[t,0]],\
#                               [e.trkdtimes[t,1]],[e.trkdtimes[t,2]],[0])
#             if len(ix)>1: print 'We have a problem'
#             elif len(ix)==0:
#                 if t==0: print 'No time match',e.trkdtimes[t]
#                 continue
#             ch = e.blobs[refkey]['ch'][e.trk[t]]
#             chmask = my.poly2mask(xypts[0],xypts[1],ch)
#             r=np.ma.MaskedArray(rain[ix,:,:],mask=~chmask)
#             rtmp[t,:,:]=r
#         rtmp=np.ma.squeeze(rtmp)
#         routput.append(rtmp)
#     routput=np.ma.asarray(routput)
#
#     return routput

def grideventmask(allpolys,s,lon,lat):
    '''Returns a mask for gridded dataset for each day of each event
    format: dictionary for each entry'''
    keylist = allpolys.keys()
    keylist.sort()

    gridmasks = {}
    masklist=[]
    maskdtimes=[]
    maskekeys=[]
    for k in keylist:
        polylist=[]
        hrsince = s.events[k].trkarrs[s.events[k].refkey][:,0]
        dtimes = my.hrnum2dates(hrsince)
        chulltrue = s.events[k].trkarrs[s.events[k].refkey][:,1] != -1
        print chulltrue, dtimes
        cnt=0
        for poly in allpolys[k]:
            if chulltrue[cnt]:
                chmask = my.poly2mask(lon,lat,poly)
                polylist.append(chmask)
                masklist.append(chmask)
                maskdtimes.append(dtimes[cnt])
                maskekeys.append(k)
            cnt+=1
        gridmasks[k] = polylist

    return gridmasks, masklist, maskdtimes, maskekeys

def ploteventrain(m,k,event,rainmasklist,raindata,griddata=False,key=False):
    '''Plot rainfall for each event based on ouput from above functions'''
    #mycm = plt.cm.GnBu
    #mycm = plt.cm.BuPu
    mycm = plt.cm.jet
    rain,dates,xypts = raindata
    #rain = np.where(rain<5, np.nan, rain)
    for it in xrange(rain.shape[0]):
        rain[it,:] = np.where(rain[it,:] < 5,np.NaN,rain[it,:])
    #olr = np.where(olr<260, np.nan, olr)
    x,y = xypts[:,0], xypts[:,1]
    ixt = []
    for t in event.trkdtimes:
        ix = my.ixdtimes(dates,[t[0]],[t[1]],[t[2]],[])
        ixt.append(ix[0])
    cnt=0
    for mask in rainmasklist:
        trkdtime = event.trkdtimes[cnt]
        print trkdtime
        plt.title('Cloudband Rainfall:' +str(trkdtime))
        if griddata:
            olat,olon,dtime,olr = griddata
            ox, oy = np.meshgrid(olon,olat)
            ixgrid = my.ixdtimes(dtime,[trkdtime[0]],[trkdtime[1]],\
                                 [trkdtime[2]],[])[0]
            m.pcolor(ox,oy,olr[ixgrid,:,:],cmap=plt.cm.PuBu_r)
            plt.clim(150,300)
        # PLOT THE CLOUD BAND CONTOUR
        plt.plot(event.blobs[event.refkey]['ch'][event.trk[cnt]][:,0],\
                 event.blobs[event.refkey]['ch'][event.trk[cnt]][:,1],\
                 color='b',ls='-',lw='3.0')
        # PLOT CONTOUR OF BLOB PROVIDED BY "KEY" WORD
        if key:
            e=event
            iarr = np.ndarray((0,),dtype=np.int32)
            hrtm = e.trktimes[cnt]
            ir = np.where((e.trkarrs[e.refkey][:,1] != -1) & \
                          (e.trkarrs[e.refkey][:,0]==hrtm))
            iarr = np.append(iarr,ir)
            # Get the blob indices for tracks we want contours for
            trkarr = np.int32(e.trkarrs[key][iarr,:])
            if trkarr.ndim==2:
                ixx = np.where(trkarr[:,1]>0)
                itrk = trkarr[ixx,1].squeeze()
            elif trkarr.ndim==3:
                itrk = np.ndarray((0,),dtype=np.int32)
                for d in xrange(trkarr.shape[2]):
                    ixx = np.where(trkarr[:,1,d]>0)
                    itrk = np.append(itrk,trkarr[ixx,1,d].squeeze())
            for ixtrk in itrk:
                plt.plot(event.blobs[key]['ch'][ixtrk][:,0],\
                         event.blobs[key]['ch'][ixtrk][:,1],\
                         color='m',ls='--',lw='2.0')

        # skip plotting if there is no data
        if ~ np.any(mask) or ~ np.any(~ np.isnan(rain[mask,ixt[cnt]])):    
            plt.text(0.5,0.5,"No rainfall inside contour")
            continue
        m.scatter(x[mask],y[mask],5,rain[mask,ixt[cnt]],\
                 cmap=mycm,edgecolor='none');plt.clim(5,100)
        m.scatter(x[~mask],y[~mask],5,np.where(rain[~mask,ixt[cnt]]<5,\
                 np.nan,rain[~mask,ixt[cnt]]),\
                 cmap=plt.cm.GnBu,edgecolor='none');plt.clim(5,50)
        ax = plt.gca()
        #if cnt==0:
            #plt.colorbar()
        m.drawcoastlines()
        m.drawcountries()
        cnt+=1
        dstring='blobs'+str(event.trkkey)+'-%d%02d%02d'\
                %(trkdtime[0],trkdtime[1],trkdtime[2])
        fname = rainfigdir+dstring+'_ee.png'
        plt.savefig(fname,dpi=200)
        #time.sleep(0.2)
        plt.axes(ax)
        plt.cla()
    #raw_input()
    #plt.clf()

def eventrainbystation(rain,dates,masklist,maskdtimes):
    erain=rain.copy()
    ixlist=[]#;edates=[]
    for i in xrange(len(maskdtimes)):
        dt = maskdtimes[i];rainmask = masklist[i]
        ix = my.ixdtimes(dates,[dt[0]],[dt[1]],[dt[2]],[])
        if len(ix)>1:print ix
        ixlist.append(ix)#;edates.append(dt)
        erain[:,ix] = np.where(rainmask[:,np.newaxis],erain[:,ix],np.nan)
    #edates=np.asarray(edates)
    return erain

def seasoncbcountmasks(s,eventkeys,yrs,lat,lon,key='noaa-olr-0-all',\
    season='coreseason',flagonly=False):
    ''' Use spatiofreq2 to return cloud band grid point counts
    USAGE: yrs is list of starting year of seasons of interest'''
    dset, varstr, levsel, deriv, exp = key.split('-')
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season
    mst,mend=mns[0],mns[-1]
    dst=1;dend=monthends[mend-1]
    if dset=='hadam3p':dend=30
    if dset=='hadam3p' and mend==2:dend=28
    cbcount=np.ndarray((0,len(lat),len(lon)),dtype=np.int32)
    for yr in yrs:
        seaskeys=specificseasons(s,eventkeys,[yr],startd=[mst,dst],\
                                                  endd=[mend,dend])
        gridcount = spatiofreq2(False,s,lat,lon,[yr],seaskeys,key=key,\
                                flagonly=flagonly)
        cbcount=np.append(cbcount,gridcount[np.newaxis,:,:],axis=0)

    return cbcount, yrs
