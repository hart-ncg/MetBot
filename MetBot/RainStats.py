'''RainStats.py: A module of the MetBot package

This module contains a collection of functions to calculate various statistics
which describe contribution of MetBot systems to given rainfall data set

Much of this is based on a masked array of event rain:
        time x stations (mask=True outside basetrack blobs' contours)
Station rainfall is included but is of the format stations x time'''
import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import mytools as my
import mynetcdf as mync
import MetBlobs as blb
import SynopticAnatomy as sy
import SynopticPlot as syp
import EventStats as stats
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from time import sleep as sleep
typlt=type(plt)
tymp=bm.Basemap

#figdir='/home/neil/work/computervision/metblobs/SA/statsfigs/'
figdir='/home/neil/work/computervision/metblobs/SA/statsfigs/rainstats/'
season=[8,9,10,11,12,1,2,3,4,5,6,7]
coreseason=[10,11,12,1,2,3,]
monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul']
mndict={}
for i in xrange(len(season)): mndict[season[i]]=monthstr[i]

def projectrainmask_in2_rawdata(rainmasks,rawrain,rawhrtime,flagdetail,everywhere=False):
    '''projectrainmask_in2_rawdata(rainmask,rawrain)

    Usage: rainmask would be output from SynopticPlot.multirainmasks
           rawrain is for example WRC time x stations data'''
    erain,etime,ekeys,flagarr = rainmasks['rain']
    # SUBSET DATA FURTHER IF SPECIFY FLAGDETAILS
    ixsub = syp.flagdetailsub(flagarr,ekeys,flagdetail)
    if len(ixsub)>0:
        erain,etime,ekeys,flagarr = erain[ixsub,:],etime[ixsub],ekeys[ixsub],flagarr[ixsub,:]
    # CORE OF THE FUNCTION
    erain = erain.T
    projdata = rawrain.copy()
    projmsk = np.ones(projdata.shape,dtype=np.bool)
    for i in xrange(rawhrtime.shape[0]):
        ix =  np.where(etime==rawhrtime[i])[0]
        if len(ix)==0:      # NO TIME MATCHES
            continue
        elif len(ix)>4:     # NUMEROUS TIME MATCHES IN MULTIPLE EVENTTRACKS
            iux, ix_iux = np.unique(ekeys[ix],return_index=True)
            #print ekeys[ix[ix_iux]]
            masklist=[]
            for index in ix_iux:
                if everywhere and np.any(~erain.mask[:,ix[index]]): erain.mask[:,ix[index]][:]=False
                #print ekeys[ix[index]]
                projmsk[:,i] = ~( (~projmsk[:,i]) | (~erain.mask[:,ix[index]]) )
        else:              # STANDARD EXPECTED FOR 4 TIME MATCHES (REASON IS HAVE 4 TIME DAILY TIMESTEPPING DUE TO NCEP SO FILL UP THESE STEPS FOR DAILY VARIABLES SUCH AS RAINFALL)
            ixt=ix[0]
            if everywhere and np.any(~erain.mask[:,ixt]): erain.mask[:,ixt][:]=False
            projmsk[:,i] = erain.mask[:,ixt]

    projected = np.ma.MaskedArray(projdata,mask=projmsk,dtype=np.float32)

    return projected

def clim_wetdayspercent(erainmasked,rtime,wetness=1):
    '''clim_wetdayspercent(erainmasked,rtime)

    Usage: erainm is outputfrom projectrainmask_in2_rawdata
           rtime is associated time
    Returns: wetdayscount (stations x 13) (ndarray), 12 months and annual mean
             freq_wet      frequency of wewet
             scycle_wet & scycle_freqwet (stations x nseasons x nmonths),
             which are same as above stats but per season'''
    season=[8,9,10,11,12,1,2,3,4,5,6,7]
    dates = my.hrnum2dates(rtime)
    years = np.unique(dates[:,0])

    wetdayscnt = np.zeros((erainmasked.shape[0],len(season)))
    freqwet = wetdayscnt.copy()
    eventmask = ~erainmasked.mask
    wetmask = erainmasked.data > wetness

    # LOOK AT THE TOTAL PERCENTAGE IN ALL DATA
    cnt=0
    for mn in season:
        ix = np.where(dates[:,1]==mn)[0]
        wetmsk = np.int32(wetmask[:,ix])
        ewetmsk = np.int32(wetmask[:,ix] & eventmask[:,ix])
        wetdayscnt[:,cnt] = wetmsk.sum(1)/np.float32(len(years))
        freqwet[:,cnt] = ewetmsk.sum(1)/np.float32(len(years))/np.float32(wetdayscnt[:,cnt])
        cnt+=1
    # SPLIT UP BY ACTUAL SEASONS
    scycle_wet = np.zeros((erainmasked.shape[0],len(years[:-1]),len(season)))
    scycle_freqwet = scycle_wet.copy()
    for iyr in xrange(len(years[:-1])):
        yr=years[iyr]
        for imn in xrange(len(season)):
            mn=season[imn]
            if imn < 5:
                ix = np.where((dates[:,0]==yr) & (dates[:,1]==mn))[0]
                wetmsk = np.int32(wetmask[:,ix])
                ewetmsk = np.int32(wetmask[:,ix] & eventmask[:,ix])
            else:
                ix = np.where((dates[:,0]==yr+1) & (dates[:,1]==mn))[0]
                wetmsk = np.int32(wetmask[:,ix])
                ewetmsk = np.int32(wetmask[:,ix] & eventmask[:,ix])
            scycle_wet[:,iyr,imn] = wetmsk.sum(1)
            scycle_freqwet[:,iyr,imn] = ewetmsk.sum(1)/np.float32(wetmsk.sum(1))
    return wetdayscnt, freqwet, scycle_wet, scycle_freqwet

def quicksubplot(m,rain,rlat,rlon,reftime,maskval=0.0,clim=[0,50,5],wait=1,saveto=False):
    plt.ion()
    plt.show()
    mycm=my.goodraincm()
    #mycm=plt.cm.gray_r
    mycm=plt.cm.jet
    if isinstance(rain,np.ma.MaskedArray):
        msk = (rain.data<=maskval) | erainm.mask
        rainm = np.ma.MaskedArray(rain.data,mask=msk,dtype=np.float32)
    else:
        msk = (rain<=maskval) | np.isnan(rain)
        rainm = np.ma.MaskedArray(rain,mask=msk,dtype=np.float32)
    #m, f = blb.SAfrBasemap(np.arange(-38,-18),np.arange(10,40),drawstuff=True,prj='cyl',fno=1,rsltn='l')
    #m, f = blb.SAfrBasemap(rlat,rlon,drawstuff=True,prj='cyl',fno=1,rsltn='l')
    for i in xrange(rainm.shape[1]):
        if reftime=='season': ttl=mndict[season[i]];exec('plt.subplot(3,4,'+str(i+1)+')')
        elif reftime=='coreseason': ttl=mndict[coreseason[i]];exec('plt.subplot(2,3,'+str(i+1)+')')
        else: ttl=str(reftime[i]);exec('plt.subplot(3,4,'+str(i+1)+')')

        rdata=rainm[:,i]
        if not np.any(~rdata.mask):
            print 'NOTHING TO SHOW'
            #m.scatter(rlon,rlat,15,rainm[:,i],cmap=plt.cm.OrRd,edgecolor='none');plt.clim(0,50)
            syp.redrawmap(m)#,lns=True,resol='hi')
            plt.text(22,-28,'NOTHING TO VALID SHOW',fontsize=5,rotation=25)
            plt.title(ttl)
            plt.draw()
            sleep(wait)
            continue
        #rdata=np.where(rdata<1,np.nan,rdata)
        #m.scatter(rlon,rlat,15,rainw[:,i],cmap=plt.cm.OrRd,edgecolor='none');plt.clim(0,50)#;plt.colorbar()
        #print np.any(~np.isnan(rdata.data))
        m.scatter(rlon,rlat,5,rdata,cmap=mycm,edgecolor='none');plt.clim(clim[:2])
        syp.redrawmap(m)#,lns=True)
        plt.title(ttl)
        plt.subplots_adjust(left=0.05,right=0.88,top=0.95,bottom=0.02,wspace=0.02,hspace=0.1)
        if i==3:
            #f=plt.gcf()
            im_col=plt.gci()
            #bounds = np.arange(clim[0],clim[1],1)
            #axcol = f.add_axes([0.90,0.2,0.015,0.5])
            #plt.colorbar(mappable=im_col,cax=axcol,boundaries=bounds)
    f=plt.gcf()
    bounds = np.arange(clim[0],clim[1]+clim[2],clim[2])
    axcol = f.add_axes([0.90,0.15,0.015,0.6])
    plt.colorbar(mappable=im_col,cax=axcol,boundaries=bounds)
    plt.draw()
    sleep(wait)
    #plt.clf()
    if saveto:
        fname=figdir+saveto+'.png'
        plt.savefig(fname,dpi=150)

def quickplot(m,rain,rlat,rlon,reftime,maskval=0.0,clim=[0,50],wait=1,saveto=False):
    plt.ion()
    plt.show()
    if isinstance(rain,np.ma.MaskedArray):
        msk = (rain.data<=maskval) | erainm.mask
        rainm = np.ma.MaskedArray(rain.data,mask=msk,dtype=np.float32)
    else:
        msk = (rain<=maskval) | np.isnan(rain)
        rainm = np.ma.MaskedArray(rain,mask=msk,dtype=np.float32)
    #m, f = blb.SAfrBasemap(np.arange(-38,-18),np.arange(10,40),drawstuff=True,prj='cyl',fno=1,rsltn='l')
    #m, f = blb.SAfrBasemap(rlat,rlon,drawstuff=True,prj='cyl',fno=1,rsltn='l')
    for i in xrange(rainm.shape[1]):
        rdata=rainm[:,i]
        if reftime=='season': ttl=mndict[season[i]]
        else: ttl=str(reftime[i])
        if not np.any(~rdata.mask):
            print 'NOTHING TO SHOW'
            #m.scatter(rlon,rlat,15,rainm[:,i],cmap=plt.cm.OrRd,edgecolor='none');plt.clim(0,50)
            syp.redrawmap(m,lns=True,resol='hi')
            plt.text(22,-28,'NOTHING TO VALID SHOW',fontsize=20,rotation=25)
            plt.title(ttl)
            plt.draw()
            sleep(wait)
            plt.clf()
            continue
        #rdata=np.where(rdata<1,np.nan,rdata)
        #m.scatter(rlon,rlat,15,rainw[:,i],cmap=plt.cm.OrRd,edgecolor='none');plt.clim(0,50)#;plt.colorbar()
        #print np.any(~np.isnan(rdata.data))
        m.scatter(rlon,rlat,5,rdata,cmap=my.goodraincm(),edgecolor='none');plt.clim(clim);plt.colorbar()
        syp.redrawmap(m,lns=True,resol='hi')
        plt.title(ttl)
        plt.draw()
        sleep(wait)
        plt.clf()
        if saveto:
            fname=figdir+saveto+'_'+ttl+'.png'
            plt.savefig(fname,dpi=150)

def eventanddaycount(s,rainmasks):
    '''Simple count of number of events and rain days per month'''
    erain,etime,ekeys,flagarr = rainmasks['rain']
    uky, ixk = np.unique(ekeys,return_index=True)
    erain,etime,ekeys,flagarr=erain[ixk,:],etime[ixk],ekeys[ixk],flagarr[ixk,:]
    ixwet = np.any(~erain.mask).nonzero()
    edates=my.hrnum2dates(etime);yrs=np.unique(edates[:,0]);yrs=yrs[:-1]

    eventcnt=np.zeros((len(yrs),len(season)),dtype=np.float32)
    daycnt=np.zeros((len(yrs),len(season)),dtype=np.float32)
    iyr=-1
    for yr in yrs:
        iyr+=1
        imn=-1
        for mn in season:
            imn+=1
            ix = np.where((edates[:,0]==yr) & (edates[:,1]==mn))
            if len(ix)>0: ix=ix[0]
            else: eventcnt[iyr,imn], daycnt[iyr,imn] = 0, 0;continue
            mnkeys = ekeys[ix];nky = len(mnkeys)
            dcnt=0
            for k in mnkeys: dcnt+=len((~np.isnan(s.events[k].rainfall['wrc'][:,0])).nonzero()[0])
            eventcnt[iyr,imn]=nky
            daycnt[iyr,imn]=dcnt

    octmarecnt=eventcnt.sum(1)
    octmardcnt=daycnt.sum(1)
    summarystats = {'eventcnt_mean': eventcnt.mean(0),'eventcnt_std': eventcnt.std(0),
               'daycnt_mean': daycnt.mean(0),'daycnt_std': daycnt.std(0),
               'ONDJFevent_mean':octmarecnt.mean(),'ONDJFevent_std':octmarecnt.std(),
               'ONDJFday_mean':octmardcnt.mean(),'ONDJFday_std':octmardcnt.std()}

    return eventcnt, daycnt, octmarecnt,octmardcnt, summarystats

def extremerainfall(s,kchoice,rain,rtime,perc=.95,areacrit=.1,selectmonths=False):
    '''extremeindices, assocTTT,assocTCOL = extremerainfall(s,kchoice,rain,rtime,perc=.95,areacriteria=.1,selectmonths=False))

    Big and ungainly approach to get out extreme rainfall day indices and associations to TTTs
    Usage: pretty obvious
    Returns: extremeindices, ndarray of extreme days in rain data set
             assocTTTs, ndarray of (extremeeventkey, ttteventkey, eventsoverlap_in_days, cutofflow_or_not)
             assocTCOL, same as above but where cutofflow_or_not is true'''
    rsort=rain.copy()
    rsort.sort(axis=1)
    ### CALCULATE EXTREME VALUE CRITERIA FOR EACH STATION
    wet = (rsort > 1.)
    extrmvals=np.zeros((rsort.shape[0]))
    for istn in xrange(rsort.shape[0]):
        iw=wet[istn,:].nonzero()[0]
        iex=int(len(iw)*perc)
        iexwet=iw[iex]
        extrmvals[istn]=rsort[istn,iexwet]

    ### PRODUCE TIME INDICES OF EXTREME DATES AND A EXTREME ARRAY MASK OF DATES
    extrm_ixt=[]
    exmsk = np.ones(rain.shape,dtype=np.bool)
    nstns=rain.shape[0]
    areacriteria=int(nstns*areacrit)
    for ixt in xrange(rain.shape[1]):
        exmsk[:,ixt] = (rain[:,ixt]>=extrmvals)
        if len(exmsk[:,ixt].nonzero()[0])>=areacriteria:
            extrm_ixt.append(ixt)
    ### GET EXTREME RAINFALL DATES
    extrm_ixt=np.asarray(extrm_ixt)
    ### IF ONLY SELECTING EXTREMES FROM PARTICULAR MONTHS
    if selectmonths:
        ixselect=[]
        exdates=my.hrnum2dates(rtime[extrm_ixt])
        for mn in selectmonths:
            ixselmn = np.where(exdates[:,1]==mn)
            if len(ixselmn)>0:ixselect.extend(ixselmn[0])
            else:continue
        extrm_ixt = extrm_ixt[ixselect]
    exrtime=rtime[extrm_ixt]

    consecutive=np.array([0],dtype=np.int32)
    consecutive=np.append(consecutive,extrm_ixt[1:]-extrm_ixt[:-1])

    extrmevents={}
    ix=0
    while ix<len(extrm_ixt):
        ixlist=[extrm_ixt[ix]]
        if ix+1>=len(extrm_ixt):extrmevents[int(rtime[ixlist[0]])]=ixlist;ix+=1;continue
        diff=consecutive[ix+1]
        if diff==1:
            ix+=1
            while diff==1:
                ixlist.append(extrm_ixt[ix])
                ix+=1
                if ix>=len(extrm_ixt):break
                diff=consecutive[ix]
        else:
            ix+=1
        extrmevents[int(rtime[ixlist[0]])]=ixlist


    ttttime,tttkeys=[],[]
    for k in kchoice:
        ktimes=list(s.events[k].trktimes)
        ttttime.extend(list(ktimes))
        for dum in xrange(len(ktimes)):tttkeys.append(k)

    ttttime,tttkeys=np.asarray(ttttime),np.asarray(tttkeys)
    extrmkeys=extrmevents.keys()
    extrmkeys.sort()
    assoc=[]
    for k in extrmkeys:
        tttk=[]
        exhrs=rtime[extrmevents[k]]
        #print exhrs
        for hr in exhrs:
            imatch = np.where(ttttime==hr)[0]
            #print int(hr),imatch,len(imatch)
            if len(imatch)==0:
                continue
            else:
                imatch=list(imatch)
                #print imatch
                tttk.extend(tttkeys[imatch])
        tkey = list(np.unique(np.asarray(tttk)))
        col=0
        if 'COL' in s.events[kchoice[0]].trkarrs.keys():
            for tk in tkey:
                if len(s.events[tk].trkarrs['COL'])>0:col+=1
        #print tkey
        ### USE THIS SETUP TO CAPTURE THE SOMETIMES 2 OR MORE CLOUD BAND EVENTS ASSOCIATED TO EXTREME RAINFALL
        #if len(tttk)>0:assoc.append((k,tkey,len(tttk)))
        #else:assoc.append((k,[0],0))
        if len(tttk)>0:assoc.append((k,tkey[-1],len(tttk),col))
        else:assoc.append((k,0,0,0))
    assoc=np.asarray(assoc,dtype=np.int32)

    assocTTT=assoc[assoc[:,1].nonzero()[0],:]
    assocTCOL=assoc[assoc[:,3].nonzero()[0],:]
    print assoc.shape,assocTTT.shape,assocTCOL.shape

    return extrm_ixt, extrmevents, assocTTT, assocTCOL

def rainprops(s,kchoice,rp,propall,rdset='wrc',pdf_title=False):
    '''propvals = rainprops(s,kchoice,rp,pdf=False)

    Return the choosen summary rainfall property for all event keys selected by kchoice
    rp (rainfall property statistic) are two parts val_stat
       where val can be: rmn, rmx, wet, hev
       and      stat is: all, mn, mx

       can plot result as a pdf'''

    rpx={'rmn_all':[5,15,1.],'rmn_mx':[5,15,1.],'rmn_mn':[5,15,1.],
         'rmx_all':[50,200,20],'rmx_mx':[100,200,10],'rmx_mn':[50,150,10],
         'wet_all':[.1,.9,.075],'wet_mx':[.1,.9,.075],'wet_mn':[.1,.9,.075],
         'hev_all':[0.,.4,.05],'hev_mx':[0.,.4,.05],'hev_mn':[0.,.4,.05],
         'lngth_lngth':[1,10,1]}
    rpy={'rmn':[0,.25,.05],'rmx':[0,.025,.005],'wet':[0,5,1],'hev':[0,9,2],'lngth':[0,.35,.05]}
    rcode={'rmn':0,'rmx':1,'wet':2,'hev':3,'lngth':0}

    rpval,rps = rp.split('_')
    rprop=[]
    for k in kchoice:
        e=s.events[k]
        eventrain=e.rainfall[rdset]
        # Check event had rainfall in dataset
        israin, notnan = len(eventrain)>0, np.any(~np.isnan(eventrain))
        if israin and notnan:
            erain=eventrain[:,rcode[rpval]]
            ixr=np.where(~np.isnan(erain))[0]
            erain=erain[ixr]
            if rps=='all':data=list(erain);rprop.extend(data)
            elif rps=='mx':data=np.nanmax(erain);rprop.append(data)
            elif rps=='mn':data=sps.nanmean(erain);rprop.append(data)
            elif rps=='lngth':data=len(ixr);rprop.append(data)
        else: rprop.append(np.nan)
    rprop=np.asarray(rprop)

    if pdf_title:
        inonan=np.where(~np.isnan(rprop))[0];rprop=rprop[inonan]
        xl=rpx[rp];yl=rpy[rpval]
        pp,xx,ss=plt.hist(propall,bins=np.arange(xl[0],xl[1],xl[2]),normed=True,align='mid',facecolor="0.3")
        plt.cla();xx=(xx[:-1]+xx[1:])/2.0
        plt.hist(rprop,bins=np.arange(xl[0],xl[1],xl[2]),normed=True,align='mid',facecolor="0.6")
        plt.plot(xx,pp,'k--',lw=2.)
        plt.xlim(xl[:2])
        plt.ylim(yl[:2]);plt.yticks(np.arange(yl[0],yl[1],yl[2]))
        plt.title(pdf_title+' (n='+str(len(rprop))+')',fontweight='normal')
        plt.savefig(pdf_title + '.' + rp + '.histogram.test.png', dpi=150)

    return rprop
