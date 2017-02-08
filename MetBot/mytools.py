# mytools.py
''' A module containing some time saving functions I use often'''

import numpy as np
try:
    import Scientific.IO.NetCDF as kh
    from mpl_toolkits.basemap import date2index, date2num, num2date
    print "Importing Scientific.IO.NetCDF library for netcdf."
    print "mpl_toolkits.basemap used for nc date support."
except ImportError:
    print "Importing netCDF4 library for netcdf & date support."
    import netCDF4 as kh
    kh.NetCDFFile = kh.Dataset
    date2index,date2num,num2date = kh.date2index,kh.date2num,kh.num2date
import scipy.stats as sps
from scipy.signal import convolve
from scipy import mgrid
import time as tm
import mynetcdf as mync
from pyclimate import diffoperators as diff
#import matplotlib.nxutils as nx 
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter
import scipy.ndimage as nim
from scipy.io import loadmat
import glob, os, socket
import cPickle
import datetime
import matplotlib.colors as colors
from matplotlib.collections import LineCollection
import errno

monthends = [31,28,31,30,31,30,31,31,30,31,30,31]
monthends_leap = [31,29,31,30,31,30,31,31,30,31,30,31]
season=[8,9,10,11,12,1,2,3,4,5,6,7]
coreseason=[10,11,12,1,2,3]
monthstr=['Aug','Sept','Oct','Nov','Dec','Jan','Feb','Mar','Apr',\
'May','Jun','Jul']
mndict={}
for i in xrange(len(season)): mndict[season[i]]=monthstr[i]

# DICTIONARY IMPORTS
from climate_index_definitions import * # gives indexdict
def country_polygon(country):
    '''lonlat_poly = country_polygon(country)
    Lifted out of matlab mapping toolbox
    worldhi.mat files of Political boundaries
    '''
    polbounds = '/home/neil/data/Topo/worldhi.mat'
    wh = loadmat(polbounds)
    C = wh[country][0][0]
    long,lat = C.long, C.lat
    ix=np.where(np.isnan(long))[0][0]
    if country=='SouthAfrica':
        long,lat = long[:ix],lat[:ix]
    return np.hstack((long,lat))

def runningmean(data,pts=3,dropVals=False):
    '''Running mean of data
    '''
    if pts%2 ==0:
        print 'pts should be an odd number!'
        return 7
    smoothed=[]
    ixoff=(pts-1)/2
    for i in xrange(len(data)-(pts-1)):
        ix=i+ixoff
        spoint=data[ix-ixoff:ix+ixoff+1].mean()
        smoothed.append(spoint)
    rmn = np.asarray(smoothed)
    if dropVals: return rmn
    if pts==3: smooth = np.hstack((rmn[0],rmn,rmn[-1]))
    if pts==5: smooth = np.hstack((rmn[0],rmn[0],rmn,rmn[-1],rmn[-1]))
    if pts==7: smooth = np.hstack((rmn[0],rmn[0],rmn[0],rmn,\
                              rmn[-1],rmn[-1],rmn[-1]))
    return smooth

def smoothgauss(vrb,smth=3):
    '''Smooths data assuming time is first dimension in a time,,lat,lon file'''
    print 'Using gaussian smoothing with std=',smth
    for i in xrange(vrb.shape[0]):
        vrb[i,:,:]=nim.gaussian_filter(vrb[i,:,:],smth)
    return vrb

def euclid(v1,v2):
    '''Squared Euclidean Distance between two scalars or equally matched vectors
       USAGE: d = euclid(v1,v2)
    '''

    d2 = np.sqrt(np.sum((v1-v2)**2))
    return d2

def AFclasses(data,cls):
    '''Function that takes data and output for K. Hinsen's ScientificPython
       affinity propagation algorithm and assigns each observation its classes
    
    Usage: classes = AFclasses(data,cls)
    Input: data (array) obs x points, raw data that was passed to
                      Scientific.Clustering.AfffinityPropagation.DataSet object
           cls (list of arrays) output of DataSet.findclusters
    Returns: classes (arr) obs x 1 contain the class of each observation
    '''

    nclasses = len(cls)
    classes = np.zeros((data.shape[0],1))
    for i in xrange(nclasses):
        iclass = cls[i]
        for j in xrange(len(iclass)):
            jx=np.where(data==iclass[j])[0][0]
            classes[jx] = i

    return classes

def clustermeans(data,clusters,obs_time):
    '''Calculates means of clusters

    Usage: centroids, cobs_time = clustermeans(data,clusters,obs_time)
    '''
    nclasses = len(np.unique(clusters))
    centr = np.zeros((nclasses,data.shape[1],data.shape[2]))
    cobs_time = []
    for i in xrange(nclasses):
        ix = np.where(clusters == i)[0]
        centr[i,:,:] = np.mean(data[ix,:,:],0)
        cobs_time.append(obs_time[ix])

    return centr, cobs_time

def regress_residual(index,field):
    '''resid_field, regressionslope = regress_residual(index,field)

    Linear Regression of each gridpoint through time with the given index to
    "remove influence" of given index.
    Returns the residual anomalies

    Input: index (array) obs x 1
                             field (array) obs x dim x dim
    Returns: resid_field (array) matches input field dimensions
    '''
    nobs,dim1,dim2 = field.shape
    expected_field = np.zeros((nobs,dim1,dim2))
    regressionslope = np.zeros((dim1,dim2))
    for i in xrange(dim1):
        for j in xrange(dim2):
            regress = sps.linregress(field[:,i,j],index)
            expected_field[:,i,j]=regress[1]+regress[0]*field[:,i,j]
            regressionslope[i,j] = regress[0]
    return field-expected_field, regressionslope


def climate_index(timemonth,data,lat,lon,ci_name):
    '''clim_index = climate_index(timemonth,data,lat,lon,ci_name)

    A simple function that returns standardized climate index.
    Contains most main indices, see climate_index_defintions.py for 
    dictionary of all supported indices

    Usage: timemonth is a list array of months of length of data
           data (time x lat x lon)
           ci_name must be a string

    indices: nino12, nino3, nino34, nino4'''

    lt1,lt2,ln1,ln2 = indexdict[ci_name]
    ilats = np.where((lt1<=lat) & (lat<=lt2))[0]
    ilons = np.where((ln1<=lon) & (lon<=ln2))[0]
    index = sps.nanmean(data[:,ilats[0]:ilats[-1]+1,ilons[0]:ilons[-1]+1],2)
    index = sps.nanmean(index[:],1)
    for i in xrange(12):
        ix = np.where(timemonth[:] == i+1)
        index[ix] = sps.zscores(index[ix])

    return index

def mon_anom_field(timemonth,data):
    '''data_anoms = mon_anom_field(timemonth,data)

             Returns the anomalies of the field, provided data has observations as first dimension'''

    t,x,y = data.shape
    data_mean = np.zeros((12,x,y))
    data_std = np.zeros((12,x,y))
    for i in xrange(12):
        ix = np.where(timemonth[:] == i+1)
        data_mean[i,:,:] = np.mean(data[ix[0],:,:],0)
        data_std[i,:,:] = np.std(data[ix[0],:,:],0)

    datameanarr=np.zeros(data.shape)
    datastdarr =np.zeros(data.shape)
    for i in xrange(12):
        datameanarr[np.where(timemonth[:]==i+1)[0],:,:] = data_mean[i,:,:]
        datastdarr[np.where(timemonth[:]==i+1)[0],:,:] = data_std[i,:,:]

    data_anoms = (data - datameanarr)/datastdarr

    return data_anoms

def daily_anom_field(dtime,data,*args):
    '''data_anoms = 6hourly_anom_field(dtime,data,[data_mean,data_std])

    If neither data_mean (nor data_std) are provided
    return tuple is (data_anoms, data_mean, data_std)

    If data_mean (and data_std) is provided in args the mean is not calculated
    but the anomalies (standardized) are returned.
    return tuple is (data_anoms)'''

    if dtime.ndim==1:
        dtime = hrnum2dates(dtime)
        print " I have run hrnum2dates conversion..."

    endlist=monthends
    if len(args)==0:
        print 'Calculating daily mean seasonal cycle field...'
        t,x,y = data.shape
        data_mean = np.zeros((365,x,y),dtype=np.float32)
        data_std = np.zeros((365,x,y),dtype=np.float32)
        data_anoms = data.copy()
        cnt=0
        for m in xrange(12):
            for d in xrange(endlist[m]):
                ix = ixdtimes(dtime,False,[m+1],[d+1],False)
                data_mean[cnt,:,:] = np.mean(data[ix,:,:],0)
                data_std[cnt,:,:]       = np.std(data[ix,:,:],0)
                data_anoms[ix,:,:] = data[ix,:,:]-\
                       np.tile(data_mean[cnt,:,:][np.newaxis,:,:],(len(ix),1,1))
                cnt+=1
        # HACK FOR LEAP YEAR'S 29TH'S
        ixlp = ixdtimes(dtime,False,[2],[29],False)
        data_anoms[ixlp,:,:] = data_anoms[ixlp,:,:] -\
                     np.tile(data_mean[57,:,:][np.newaxis,:,:],(len(ixlp),1,1))
        out=(data_anoms, data_mean, data_std)
    elif len(args)>0:
        try: data_mean, data_std = args
        except: data_mean, data_std = args[0],False
        data_anoms = data.copy()
        cnt=0
        for m in xrange(12):
            for d in xrange(endlist[m]):
                ix = ixdtimes(dtime,False,[m+1],[d+1],False)
                data_anoms[ix,:,:] = data[ix,:,:]-\
                       np.tile(data_mean[cnt,:,:][np.newaxis,:,:],(len(ix),1,1))
                if isinstance(data_std,np.ndarray):
                    data_anoms[ix,:,:] = (data[ix,:,:]-\
                    np.tile(data_mean[cnt,:,:][np.newaxis,:,:],(len(ix),1,1)))/\
                    np.tile(data_std[cnt,:,:][np.newaxis,:,:],(len(ix),1,1))
                cnt+=1
    # HACK FOR LEAP YEAR'S 29TH'S
        ixlp = ixdtimes(dtime,False,[2],[29],False)
        data_anoms[ixlp,:,:] = data_anoms[ixlp,:,:] -\
          np.tile(data_mean[57,:,:][np.newaxis,:,:],(len(ixlp),1,1))
        if isinstance(data_std,np.ndarray):
            data_anoms[ixlp,:,:] = (data_anoms[ixlp,:,:] -\
               np.tile(data_mean[57,:,:][np.newaxis,:,:],(len(ixlp),1,1)))/\
               np.tile(data_std[57,:,:][np.newaxis,:,:],(len(ixlp),1,1))
        out=(data_anoms)

    return out

def sixhourly_anom_field(dtime,data,*args):
    '''data_anoms = sixhourly_anom_field(dtime,data,[data_mean,data_std])

    If neither data_mean (nor data_std) are provided
    return tuple is (data_anoms, data_mean, data_std)

    If data_mean (and data_std) is provided in args the mean is not calculated
    but the anomalies (standardized) are returned.
    return tuple is (data_anoms)'''

    if dtime.ndim==1:
        dtime = hrnum2dates(dtime)
        print " I have run hrnum2dates conversion..."

    endlist=[31,28,31,30,31,30,31,31,30,31,30,31]
    if len(args)==0:
        print 'Calculating sixhourly mean seasonal cycle field...'
        t,x,y = data.shape
        data_mean = np.zeros((1460,x,y),dtype=np.float32)
        data_std = np.zeros((1460,x,y),dtype=np.float32)
        data_anoms = data.copy()
        cnt=0
        for m in xrange(12):
            for d in xrange(endlist[m]):
                for h in [0,6,12,18]:
                    ix = ixdtimes(dtime,False,[m+1],[d+1],[h])
                    data_mean[cnt,:,:] = np.mean(data[ix,:,:],0)
                    data_std[cnt,:,:]       = np.std(data[ix,:,:],0)
                    data_anoms[ix,:,:] = data[ix,:,:]-\
                    np.tile(data_mean[cnt,:,:][np.newaxis,:,:],(len(ix),1,1))
                    cnt+=1
        # HACK FOR LEAP YEAR'S 29TH'S
        cnt2=0
        for h in [0,6,12,18]:
            ixlp = ixdtimes(dtime,False,[2],[29],[h])
            data_anoms[ixlp,:,:] = data_anoms[ixlp,:,:] - \
            np.tile(data_mean[231+cnt2,:,:][np.newaxis,:,:],(len(ixlp),1,1))
            cnt2+=1
        out=(data_anoms, data_mean, data_std)

    elif len(args)>0:
        try: data_mean, data_std = args
        except: data_mean, data_std = args[0],False
        data_anoms = data.copy()
        cnt=0
        for m in xrange(12):
            for d in xrange(endlist[m]):
                for h in [0,6,12,18]:
                    ix = ixdtimes(dtime,False,[m+1],[d+1],[h])
                    data_anoms[ix,:,:] = data[ix,:,:]-\
                    np.tile(data_mean[cnt,:,:][np.newaxis,:,:],(len(ix),1,1))
                    if isinstance(data_std,np.ndarray):
                        data_anoms[ix,:,:] = (data[ix,:,:]-np.tile(\
                        data_mean[cnt,:,:][np.newaxis,:,:],(len(ix),1,1)))/\
                        np.tile(data_std[cnt,:,:][np.newaxis,:,:],(len(ix),1,1))
                    cnt+=1
        # HACK FOR LEAP YEAR'S 29TH'S
        cnt2=0
        for h in [0,6,12,18]:
            ixlp = ixdtimes(dtime,False,[2],[29],[h])
            data_anoms[ixlp,:,:] = data_anoms[ixlp,:,:] -\
            np.tile(data_mean[231+cnt2,:,:][np.newaxis,:,:],(len(ixlp),1,1))
            if isinstance(data_std,np.ndarray):
                data_anoms[ixlp,:,:] = (data_anoms[ixlp,:,:] - np.tile(\
                data_mean[231+cnt2,:,:][np.newaxis,:,:],(len(ixlp),1,1)))/\
                np.tile(data_std[231+cnt2,:,:][np.newaxis,:,:],(len(ixlp),1,1))
            cnt2+=1
        out=(data_anoms)

    return out

def closestlonlat(lon,lat,pt):
    '''Returns the ilt, iln indices for gridpoint closest to the 
       (lon,lat) values provided in pt'''
    ln,lt=pt
    londiff,latdiff = np.abs(lon-ln), np.abs(lat-lt)
    lndmin, ltdmin = londiff.min(), latdiff.min()
    iln,ilt=(londiff==lndmin).nonzero()[0][0],(latdiff==ltdmin).nonzero()[0][0]
    return (iln, ilt)

def subset(data,lat,lon,sub):
    ''' subdata, sublat, sublon = subset(data,lat,lon,sub)

    Simple function to subset data

    INPUT: lat (1-dim array)
                             lon (1-dim array)
                             data (obs/time x lats x lon array)
                             sub ((latmin,latmax),(lonmin,lonmax)) tuple
                             or
                             sub can be "Fauch08", "SH", "imkf", "playing"

    RETURNS: subdata, sublat, sublon
    '''
    domains={}
    domains['Fauch08'] = ((-40.0,-15.0),(7.5,70.0))
    domains['Fauch08_mod'] = ((-40.0,-10.0),(7.5,70.0))
    domains['SH'] = ((-90.0,0.0),(0,360))
    domains['imkf'] = ((-46.25,-10.0),(0.0,91.875))
    domains['playing'] = ((-60.0,-5.0),(0.0,70.0))
    domains['SA'] = ((-60.0,0.0),(0.0,80.0))


    if isinstance(sub,str):
        domain=domains[sub]
    elif isinstance(sub,tuple):
        domain=sub
    else:
        print "ERROR: sub needs to be a tuple or a string"

    lt1,lt2 = domain[0]; ln1,ln2 = domain[1]
    ilats = ((lt1<=lat) & (lat<=lt2)).nonzero()[0]
    ilons = ((ln1<=lon) & (lon<=ln2)).nonzero()[0]
    subdata = data[:,ilats[0]:ilats[-1]+1,ilons[0]:ilons[-1]+1]
    sublat=lat[ilats[0]:ilats[-1]+1];sublon=lon[ilons[0]:ilons[-1]+1]

    return subdata, sublat, sublon

def gauss(mu, sigma, x ):
    return (1.0/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2.0*sigma**2))

def gaussian(vec,a='sn',b=0,c=1):
    '''yx = gaussian(vec,a='sn',b=0,c=1))

    Applies gaussian function to given vector

    Input: vec (array) vector of data to give as input to gaussian
           a : amplitude, defaults to standard normal distribution
           b : mean, default = 0
           c : variance, defaults = 1
    Returns: yx
    '''
    if a == 'sn':
        a = 1/(2*pi*c**2)**0.5
    A = (vec-b)**2; B = 2*c**2
    yx = a*np.exp(-A/float(B))
    return yx

def randindgen(idxlen,n_irnd):
    '''irandom = randindgen(idxlen,n_irnd)
    Random Index Generator
    Generates a set of random integers on the interval [0,idxlen)
    Useful as index of random slices from a large array'''
    irandom = np.int16(np.round(np.random.random((n_irnd))*(idxlen-1)))
    return  irandom

def CompareAllIndex(ix):
    '''ix_permutes = CompareAllIndex(ix)
    A simple function to return an array of pairs of indices to
    perform whatever comparison need be done

    Was written to be able to get distance vectors between each point on an
    image. Provides index for all permutations, doesn't actually utilise 
    symmetry inherent in the question ie array has entries both [1,0] and [0,1]
    '''
    n=ix.shape[0]
    perms=[]
    for i in xrange(n):
        for j in xrange(n):
            perms.append([i,j])
    perms=np.asarray(perms)
    ix_permutes=np.hstack((ix[perms[:,0]][:,newaxis],ix[perms[:,1]][:,newaxis]))

    return ix_permutes

def OpenMultipleNC(nclist,varstr,dset,sub=False,levselect=False):
    '''var, time, lat, lon, lev, dtime = OpenMultipleNC(nclist,varstr,dset,\
                                                     sub=False,levselect=False)
    A hack of a function to return a single variable stored across 
    mulitple .nc files:
    USAGE:
    srcpath (str): is full path prefix to files excluding month and .nc 
                   (eg. "/home/neil/data/daily/olr/olr." ; \
                    "/home/neil/data/flavours/clim.caaij.olr.")
    mnlist (list): list of 1 or more months generally (could be years) to 
                   produce filelist prefixed by srcpath
    varstr: self-explanatory
    dset (str): ncep, cfsr, hadam3p (HadAM3p experiment): just needed in 
                determining which mynetcdf.py function to use
    sub (tuple): ((latmin,latmax),(lonmin,lonmax)), default is False, 
                                                    hence no subsetting
    levselect (number): value of level to select out of data, default is False,
                        hence assumed data is time x lat x lon
    RETURNS:
    self-evident
    '''
    lev=()
    cnt=0
    for i in nclist:
        print cnt+1,": Opening ",i
        ct=str(cnt)
        if levselect:
            exec('var'+ct+', time'+ct+', lat, lon, lev, dtime'+ct+' = mync.open'+dset+'(\"'+i+'\",\"'+varstr+'\",subs=\''+sub+'\',levsel='+str(levselect)+')')
            exec('var'+ct+'=var'+ct+'.squeeze()')
        elif sub:
            exec('var'+ct+', time'+ct+', lat, lon, dtime'+ct+'= mync.open'+dset+'(\"'+i+'\",\"'+varstr+'\",subs=\''+sub+'\')')
        else:
            exec('var'+ct+', time'+ct+', lat, lon, dtime'+ct+' = mync.open'+dset+'(\"'+i+'\",\"'+varstr+'\")')
        cnt=cnt+1
    if varstr=='precipitation':
        print "Adding new axis to trmm vars..."
        for i in xrange(cnt):
            exec('var'+str(i)+'=var'+str(i)+'[np.newaxis,:,:]')
    varlist=['var'+str(i) for i in range(cnt)]
    timelist=['time'+str(i) for i in range(cnt)]
    dtimelist=['dtime'+str(i) for i in range(cnt)]
    exec('vartup=('+', '.join(varlist)+',)')
    exec('timetup=('+', '.join(timelist)+',)')
    exec('dtimetup=('+', '.join(dtimelist)+',)')
    var=np.vstack(vartup)
    time=np.hstack(timetup)
    dtime=np.vstack(dtimetup)

    return var, time, lat, lon, lev, dtime

def ixdtimes(t,yrs,mns,ds,hs, mode='match'):
    '''indices = ixdtime(times, yrs, mns, ds, hs, mode="match"/"range")
       USAGE:
       times is ntime x 4 array of entries [yr,mn,dy,hr]
       yrs, mns, ds, hs are singleton, or non-singleton lists of values match
       mode="range" returns indices of times between two times held within
                    pairs in yrs, mns, ds, hs
    '''
    if mode=='match':
        teststr='';ytest=[];mtest=[];dtest=[];htest=[]
        if yrs:
            ytest=["(t[:,0]=="+str(y)+")" for y in yrs]
            teststr=" | ".join(ytest) + " & "
        if mns:
            mtest=["(t[:,1]=="+str(m)+")" for m in mns]
            teststr=teststr + " | ".join(mtest) + " & "
        if ds:
            dtest=["(t[:,2]=="+str(d)+")" for d in ds]
            teststr=teststr + " | ".join(dtest) + " & "
        if hs:
            htest=["(t[:,3]=="+str(h)+")" for h in hs]
            teststr=teststr + " | ".join(htest)
        if teststr[-2:-1] == "&":
            teststr=teststr[:-2]
        #try:
        #print teststr
        exec('ix = np.where('+teststr+')[0]')
        #except IndexError:
            #print 'mytools.ixdtimes no valid time matches'
            #return False
    if mode=='range':
        teststr='';ytest=[];mtest=[];dtest=[];htest=[]
        tmin=yrs[0]*10e6 + mns[0]*10e4 + ds[0]*10e2 + hs[0]
        tmax=yrs[1]*10e6 + mns[1]*10e4 + ds[1]*10e2 + hs[1]

        tnum=t[:,0]*10e6 + t[:,1]*10e4 + t[:,2]*10e2 + t[:,3]
        ix = np.where((tnum >= tmin) & (tnum <= tmax))[0]

    return ix

def ixdatetimes(t,matchdates,mode='match',strict=True):
    '''indices = ixdatetimes(t,matchdates)'''
    ### ENSURE FORMATED TO NUMPY (>v1.7) DATETIME64 FORMAT
    if t.dtype==np.dtype('O'): dt64=np.asarray([np.datetime64(tt) for tt in t])
    elif t.dtype==np.dtype('datetime64[us]'): dt64 = t
    else: print 'datetimes not in understood datatype';return 0

    ### MODE TO MATCH DATES
    if mode=='match':
        ixlist=[]
        for mdt in matchdates:
            try:
                ixnew=np.where(dt64==np.datetime64(mdt))[0][0]
                ixlist.append(ixnew)
            except:
                print 'No match for', mdt
        ix = np.asarray(ixlist)
    ### MODE TO GET DATES WITHIN RANGE
    if mode=='range':
        tmin, tmax = np.datetime64(matchdates[0]), np.datetime64(matchdates[1])
        if tmax>dt64[-1]:
            print 'mytools.ixdatetimes: End date outside of datetimes range!'
            if strict: return False
        ix = np.where((dt64 >= tmin) & (dt64 <= tmax))[0]
        if len(ix)==0:
            print 'mytools.ixdatetimes: no dates match in this date range!'
            if strict: return False

    return ix

def validtimeindices(dtime,seasonstartdate,seasonlength="120 days",\
             notconsecutiveyears=False):
    ''' ixtvalid = validtimeindices(dtime,seasonstartdate,\
                    seasonlength="120 days",notconsecutiveyears=False)
        Returns: boolean array of same length of dtime which can be used
                 to indicate indices that are valid for given season time
                 period of interest.
    '''
    startmonth,startdate = seasonstartdate
    seaslen, lenunit = seasonlength.split()
    if lenunit=='days': tdelta = datetime.timedelta(days=float(seaslen))
    elif lenunit=='weeks': tdelta = datetime.timedelta(weeks=float(seaslen))
    elif lenunit=='months':
        print 'Converting to days, 30 days per month'
        tdelta = datetime.timedelta(days=float(seaslen)*30)
    yrslist=np.asarray([t.year for t in dtime])
    yrs=np.unique(yrslist)
    if isinstance(notconsecutiveyears,np.ndarray):yrs=notconsecutiveyears
    ixtbyear =[]
    for iyr, yr in enumerate(yrs):
        day1=datetime.datetime(yr,startmonth,startdate)
        daylast= day1 + tdelta
        ixt = ixdatetimes(dtime,[day1,daylast],mode='range')
        ixtbyear.append(ixt)
    ixtvalid=np.in1d(range(dtime.shape[0]),np.hstack(ixtbyear))
    return ixtvalid

def ixtwindow(href,hr1,twindow,*args,**kwargs):
    '''ixt = ixtwindow(href,hr,twindow,[dtime,data,time])
    Find the indices where times in hr1 fall within a given time window
    around the times in href.
   
    Typical example of use is to get indices into hr1 which will give the times
    close to a set of times in href which indicate hourtimes of events of
    interest.

    href and hr1 units="hours since DATE" hourtime values. Obviously,
    this then requires href an hr1 to have EXACTLY the same units.

    Returns indices of where hourtime falls within
    time window of another set of times
    Two modes to use:
    GET INDICES
    ixt = ixtwindow(href,hr,twindow)
    or
    DO SUBSETS
    ixt, [time_ixt, hgt_ixt, dtime_ixt] = 
                                   my.ixtwindow(reftime,time,24,time,hgt,dtime)
    OPTIONAL USE:
    To look only for times before or after the href times, can use
    ixt = ixtwindow(href,hr,twindow,twosided=False)
    In this case the sign of twindow will be used to determine whether to
    consider time window which extends to before or after href event times.
    '''
    if len(kwargs)>0:
        twosided=kwargs['twosided']
    else:
        twosided=True
    before, after = False, False
    if not twosided:
        if twindow>0: after=True
        elif twindow<0: before=True
    ixt=np.ndarray((0,),dtype=np.int16)
    for i, tref in enumerate(href):
        deltat = hr1 -tref 
        if twosided: ix = np.where(np.abs(deltat) <= twindow)[0]
        elif after:  ix = np.where((deltat<=twindow) & (deltat>0))[0]
        elif before:  ix = np.where((deltat>=twindow) & (deltat<0))[0]
        ixt = np.append(ixt,ix)
    ixt = np.unique(ixt)

    if len(args) > 0:
        print "Time subsetting..."
        sbstlist=[]
        for a in args:
            if a.ndim==1: sbstlist.append(a[ixt])
            elif a.ndim>1: sbstlist.append(a[ixt,...])
        return ixt, sbstlist
    else:
        return np.unique(ixt)

def hrnum2dates(hrssince,returntype='dtimearr',\
    units="hours since 1800-01-01 00:00:00.0",calendar='gregorian'):
    '''dates=hrnum2dates(hrssince,returntype='dtimearr'\
    units="hours since 1800-01-01 00:00:00.0",calendar='gregorian')
    
    Returns human readable dates either as dtimearr of enteries [YYYY,MM,DD,HH]
            or as array of datetime objects, if specify returntype="datetime"
    '''
    dtime = mync.num2date(hrssince,units=units,calendar=calendar)
    if returntype=='dtimearr':
        dates = mync.dtime2arr(dtime)
    else:
        dates = dtime
    return dates

def dates2hrnum(dates,\
              units="hours since 1800-01-01 00:00:00.0",calendar='gregorian'):
    '''hrtimes = dates2hrnum(dates,units="hours since 1800-01-01 00:00:00.0",
                              calendar='gregorian')
    Given an array of human readable dates return array of hourssince with
     given units.
    USAGE: dates can either be array of datetime.datetime objects
                  or and array of entries of format [YYYY, MM, DD, HH]'''
    if isinstance(dates[0],datetime.datetime):
        hrtime = mync.date2num(dates,units=units,calendar=calendar)
    else:
        datetimes=[]
        for dt in dates:
            y,m,d,h = dt
            datetimes.append(datetime.datetime(y,m,d,h))
        hrtime = mync.date2num(datetimes,units=units,calendar=calendar)
    return hrtime

## Matlab clones
def poly2mask(lon,lat,poly,invert=False):
    ''' maskgrid = poly2mask(lat,lon,poly,invert=False)
    
    Return a mask for values within the given polygon
    USAGE: 
    '''
    ln,lt = np.meshgrid(lon,lat)
    lls = ln.shape
    xypts = np.hstack((ln.ravel()[:,np.newaxis], lt.ravel()[:,np.newaxis]))
    #tes = np.reshape(xypts[:,0],lls)
    polypath=Path(poly)
    truelist=[]
    for xy in xypts:
        yesno = polypath.contains_point(xy)
        truelist.append(yesno)
    mask=np.asarray(truelist,dtype=np.bool8)
    maskgrid = np.reshape(mask, lls)
    if invert:
        maskgrid = ~maskgrid
    return maskgrid

def points_inside_poly(points,poly):
    thepoly=Path(poly)
    boolean=thepoly.contains_points(points)
    return boolean

## METEOROLOGICAL VARIABLES
def d_dx(var,time,lat,lon):
    '''d_dx*1e7 = d_dx(var,time,lat,lon)

    Calculates d(var)/dx of geophysical field
    using pyclimate.diff.HGRADIENT
    which calculate these values for spherical coordinates.

    Usage: Bit of hack but time is either time array or [False]

    Returns: d_dx*1e7 '''

    Grd=diff.HGRADIENT(lat,lon)
    if not time[0]:
        d_dx,d_dy=Grd.hgradient(var)
    else:
        d_dx=np.zeros(var.shape, dtype=np.float32)
        for t in xrange(len(time)):
            ugrd,vgrd=Grd.hgradient(var[t,:,:].squeeze())
            d_dx[t,:,:]=ugrd

    return d_dx*1e7

def del2(var,time,lat,lon):
    '''del2*1e10 = del2(var,time,lat,lon)

    Calculates del2 (Laplace Operator) of geophysical field
    using pyclimate.diff.HGRADIENT
                            pyclimate.diff.HDIVERGENCE
    which calculate these values for spherical coordinates.


    Usage: Bit of hack but time is either time array or [False]

    Returns: del2*1e10 '''

    Grd=diff.HGRADIENT(lat,lon)
    Div=diff.HDIVERGENCE(lat,lon)
    if not time[0]:
        ugrd,vgrd=Grd.hgradient(var)
        del2=Div.hdivergence(ugrd,vgrd)
    else:
        del2=np.zeros(var.shape, dtype=np.float32)
        for t in xrange(len(time)):
            ugrd,vgrd=Grd.hgradient(var[t,:,:].squeeze())
            del2[t,:,:]=Div.hdivergence(ugrd,vgrd)

    return del2*1e10

def div(u,v,time,lat,lon):
    '''Divergence of geophysical field on lat lon grid
    Usage: Bit of hack but time is either time array or [False]

    Returns: del2*1e10 '''

    Div=diff.HDIVERGENCE(lat,lon)
    if not time[0]:
        dv=Div.hdivergence(u,v)
    else:
        dv=np.zeros(u.shape, dtype=np.float32)
        for t in xrange(len(time)):
            dv[t,:,:]=Div.hdivergence(u[t,:,:],v[t,:,:])

    return dv

def magdir(u,v):
    '''returns the mag and bearing of vector describe by u & v'''
    mag=np.sqrt(u**2 + v**2)
    raddir=np.arctan(v/u)
    degdir=450-360*(raddir/(2*np.pi))
    return mag, degdir

def MagDir(u,v):
    print "Depreciated: Rather use mytools.magdir"
    mag, degdir = magdir(u,v)
    return mag, degdir

def coriolosscale(vrb,lat,lon,):
    '''Crude method to scale variables by coriolos force'''
    lat_grid=np.tile(lat[:,np.newaxis],(1,lon.shape[0]))
    lat_grid=np.tile(lat_grid,(vrb.shape[0],1,1))
    f0=2*7.292e-5*np.sin(-90.*np.pi/180.0)
    f=2*7.292e-5*np.sin(lat_grid*np.pi/180.0)
    scale=1-f/f0
    return vrb*scale

#def psat(T)
    #Ts             = 373.16                 # steam point temperature in K
    #ews     = 1013.246      # saturation pressure at steam point
    ## temperature, normal atmosphere
    #a1 =   5.02808
    #a2 = -1.3816E-7
    #a3 =   8.1328E-3
    #b1 = -7.90298
    #b2 = 11.344
    #b3 = -3.49149
    #psatdew = 10.^(b1*(Ts./T-1.)+a1*log10(Ts./T)+a2*(10.^(b2*(1-T./Ts))-1)\
    #+a3*(10.^(b3*(Ts./T-1))-1)+log10(ews))
    ## Ice
    #ei0 = 6.1071 # mbar
    #T0 = 273.16 # freezing point in K
    #psatice = 10.^(-9.09718*(T0./T-1.)-3.56654*log10(T0./T)+0.876793*(1.-T/T0)+log10(ei0));
def cmap_discretize(cmap, N):
     """Return a discrete colormap from the continuous colormap cmap.
         cmap: colormap instance, eg. cm.jet. 
         N: number of colors.
        Courtesy: scipy cookbook on Colormap Transformations
     """
 
     if type(cmap) == str:
         cmap = matplotlib.cm.get_cmap(cmap)
     colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
     rgbas = cmap(colors_i)
     ixs = np.linspace(0, 1., N+1)
     cdict = {}
     for ki,key in enumerate(('red','green','blue')):
         cdict[key]=[(ixs[i],rgbas[i-1,ki],rgbas[i,ki]) for i in xrange(N+1)]
     # Return colormap object.
     cm=matplotlib.colors.LinearSegmentedColormap(cmap.name+"_%d"%N,cdict,1024)
     return cm

def savef(fdir=os.getenv('HOME')+'/workingfigs/atemp_python'):
    fnum=plt.gcf().number
    figname="%s/currentfig%d.png" %(fdir, fnum)
    plt.savefig(figname)
    print 'Figure available: ', figname

## MY RAIN HANDLING ROUTINES
def goodraincm(*args):
    '''Registers a good colormap for rain, lightblue-to-darkblue-to-red
    returns object <matplotlib.colors.LinearSegmentedColormap>
    USAGE: don't pass anything to function call to return colormap
    if pass *anything* sensible  will return reversed version
    '''
    # one used I made in matlab
    _rainbluered = ((0.0, 1.0, 1.0),(0.0, 0.0, 1.0),(1.0, 0.0, 0.2))
    #_rainbluered = ((1.0, 0.95, 0.3),(0.0, 1.0, 0.8),(0.0, 0.3, 1.0),(0.4, 0.0, 1.0),(1.0, 0.0, 0.2))
    LUTSIZE=256
    cmapspec = _rainbluered
    cmapname = 'rainBuRd'
    cmapname_r = 'rainBuRd_r'
    revspec = list(reversed(cmapspec))
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    rainBuRd = colors.LinearSegmentedColormap.from_list(cmapname, \
                                                        cmapspec, LUTSIZE)
    rainBuRd_r = colors.LinearSegmentedColormap.from_list(cmapname_r,\
                                                          revspec, LUTSIZE)
    if len(args)==0: out = rainBuRd
    elif len(args)>0: out = rainBuRd_r
    return out

def betterraincm(*args):
    '''Registers a good colormap for rain, lightblue-to-darkblue-to-red
    returns object <matplotlib.colors.LinearSegmentedColormap>
    USAGE: don't pass anything to function call to return colormap
    if pass *anything* sensible  will return reversed version
    '''
    # one used I made in matlab
    _rainbluered = ((1.,1.,1.),(0.0, 1.0, 1.0),(0.0, 0.0, 1.0),(1.0, 0.0, 0.2))
    #_rainbluered = ((1.0, 0.95, 0.3),(0.0, 1.0, 0.8),(0.0, 0.3, 1.0),(0.4, 0.0, 1.0),(1.0, 0.0, 0.2))
    LUTSIZE=256
    cmapspec = _rainbluered
    cmapname = 'rainBuRd'
    cmapname_r = 'rainBuRd_r'
    revspec = list(reversed(cmapspec))
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    rainBuRd = colors.LinearSegmentedColormap.from_list(cmapname, \
                                                        cmapspec, LUTSIZE)
    rainBuRd_r = colors.LinearSegmentedColormap.from_list(cmapname_r,\
                                                          revspec, LUTSIZE)
    if len(args)==0: out = rainBuRd
    elif len(args)>0: out = rainBuRd_r
    return out

def goodwindcm(*args):
    '''Registers a good colormap for wind, lightblue-to-darkblue-to-red
    returns object <matplotlib.colors.LinearSegmentedColormap>
    USAGE: don't pass anything to function call to return colormap
    if pass *anything* sensible  will return reversed version'''
    _windbluered = ((1.0, 1.0, 1.0),(0.0, 0.0, 1.0),(1.0, 0.0, 0.2))
    #_windbluered = ((1.0, 0.95, 0.3),(0.0, 1.0, 0.8),(0.0, 0.3, 1.0),(0.4, 0.0, 1.0),(1.0, 0.0, 0.2))
    LUTSIZE=256
    cmapspec = _windbluered
    cmapname = 'windBuRd'
    cmapname_r = 'windBuRd_r'
    revspec = list(reversed(cmapspec))
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    windBuRd = colors.LinearSegmentedColormap.from_list(cmapname, \
                                                        cmapspec, LUTSIZE)
    windBuRd_r = colors.LinearSegmentedColormap.from_list(cmapname_r,\
                                                          revspec, LUTSIZE)
    if len(args)==0: out = windBuRd
    elif len(args)>0: out = windBuRd_r
    return out

def grey_color_scale(nlevels=10,rev=False):
    '''Defines a colour scale that it colourful, not too sickly and which has
       monotonically increasing shades of grey when printed in black and white 
       so that there is no ambiguity.
       Courtesy: Peter Clark and others'''
    nnxx=nlevels-1
    R=np.hstack((np.repeat(1.0,nnxx),np.arange(1.0,0.0,-1.0/nnxx),\
                 np.arange(0.0,0.75,0.75/nnxx),np.arange(0.75,1.0,0.25/nnxx),\
                 np.arange(1.0,0.0,-1.0/nnxx),np.repeat(0.0,nnxx) ))
    G=np.hstack((np.arange(1.0,0.5,-0.5/nnxx),np.arange(0.5,1.0,0.5/nnxx),\
                 np.arange(1.0,0.75,-0.25/nnxx),np.arange(0.75,0.0,-0.75/nnxx),\
                 np.repeat(0.0,nnxx),np.repeat(0.0,nnxx) ))
    B=np.hstack((np.repeat(1.0,nnxx),np.repeat(1.0,nnxx),\
                 np.arange(1.0,0.0,-1.0/nnxx),np.repeat(0.0,nnxx),\
                 np.arange(0.0,0.5,0.5/nnxx),np.arange(0.5,0.0,-0.5/nnxx) ))
    RGB=np.hstack((R[:,np.newaxis],G[:,np.newaxis],B[:,np.newaxis]))
    gcs = colors.ListedColormap(RGB,name='grey_color_scale')
    gcs_r =colors.ListedColormap(RGB[::-1,:]) 
    if not rev: out = gcs
    elif rev: out = gcs_r
  
    return out

def goodolrcm(*args):
    '''Registers a good colormap for wind, lightblue-to-darkblue-to-red
    returns object <matplotlib.colors.LinearSegmentedColormap>
    USAGE: don't pass anything to function call to return colormap
    if pass *anything* sensible  will return reversed version'''
    _olrbluered = ((1.0, 1.0, 1.0),(0.4, 0.4, 0.4),(0.8, 0.85, 1.0))
    #_windbluered = ((1.0, 0.95, 0.3),(0.0, 1.0, 0.8),(0.0, 0.3, 1.0),(0.4, 0.0, 1.0),(1.0, 0.0, 0.2))
    LUTSIZE=256
    cmapspec = _olrbluered
    cmapname = 'olrGrBu'
    cmapname_r = 'olrGrBu_r'
    revspec = list(reversed(cmapspec))
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    #pyplotmodule.cm.register_cmap(name=cmapname,data=cmapspec,lut=LUTSIZE)
    windBuRd = colors.LinearSegmentedColormap.from_list(cmapname, \
                                                        cmapspec, LUTSIZE)
    windBuRd_r = colors.LinearSegmentedColormap.from_list(cmapname_r,\
                                                          revspec, LUTSIZE)
    if len(args)==0: out = windBuRd
    elif len(args)>0: out = windBuRd_r
    return out

def trmmSAstations():
    '''Opens trmm and cuts out for the land and turns data into station like 
       format
    date x station with associate lat, lon pair for each station'''
    atrmmfile='/home/neil/data/trmm/trmm.1998.daily.SA.nc'
    trmm, time, lat, lon, dtime = mync.openncep2(atrmmfile,'trmm')
    llt,lln = np.meshgrid(lat,lon)
    lln, llt = lln.ravel()[:,np.newaxis], llt.ravel()[:,np.newaxis]
    xytrmm = np.hstack((lln, llt))
    sapol = country_polygon('SouthAfrica')
    landmask = points_inside_poly(xytrmm,sapol)
    ix = np.where(landmask)[0]
    stlon, stlat = lln[ix], llt[ix]

    trmmstations=np.zeros(())
    flist=glob.glob('/home/neil/data/trmm/trmm.*.daily.SA.nc');flist.sort()
    trmmstations=np.ndarray((ix.shape[0],0))
    dates=np.ndarray((0,4))

    for f in flist:
        print f
        trmm, time, lat, lon, dtime = mync.openncep2(f,'trmm')
        for t in xrange(trmm.shape[0]):
            trmmsts = trmm[t,:,:].ravel()[ix]
            trmmstations = np.append(trmmstations,trmmsts[:,np.newaxis],axis=1)
            dates = np.append(dates,dtime[t][np.newaxis,:],axis=0)

    trmmrain={'lat':stlat,'lon':stlon,'dates':np.int32(dates),\
              'rain':trmmstations}
    tfile='/home/neil/data/trmm/trmmSAstations2.pickle'
    pickf = open(tfile,'w')
    cPickle.dump(trmmrain,pickf)
    pickf.close()

def openWRCstations(datadir='/home/neil/sogehome/data/',\
                    f='wrcrain/wrcSA79_99.mat'):
    f=datadir+f
    print 'Opening ',f
    wrcrain = loadmat(f)
    rdate,rlon,rlat = wrcrain['dates'], wrcrain['lon'], wrcrain['lat']
    rain = wrcrain['rain']
    dates=[]
    for i in rdate:
        y,m,d=str(i[0])[:4],str(i[0])[4:6],str(i[0])[6:8]
        dates.append((int(y),int(m),int(d),0))
    dates=np.asarray(dates)
    hrtime = dates2hrnum(dates)
    out=(np.float32(rain), np.float32(hrtime), np.float32(rlat),\
          np.float32(rlon), np.int32(dates))
    return out 

def openTRMMstations(f='/home/neil/data/trmm/trmmSAstations2.pickle'):
    print 'Opening ',f
    tf = open(f,'r')
    train = cPickle.load(tf)
    rain = train['rain']
    rdate, rlat, rlon =  train['dates'], train['lat'], train['lon']
    tf.close()
    rdate[:,3]=0
    hrtime = dates2hrnum(rdate)
    out = (np.float32(rain), np.float32(hrtime), np.float32(rlat),\
           np.float32(rlon), np.float32(rdate))
    return out

def openSATOPO(f='/home/neil/data/Topo/etopo2_southernAfrica.nc'):
    #topodata = loadmat(f)
    ncf = kh.NetCDFFile(f,'r')
    topodata, lat, lon = ncf.variables['topo'][:], ncf.variables['lat'][:],\
                         ncf.variables['lon'][:]
    data = np.ma.MaskedArray(topodata,mask=topodata<0,dtype=np.float32)
    return data, lat, lon

def summerrain(data='wrc'):
    summerregion=np.asarray([(24,-25),(28,-34),(35,-33),(35,-20),\
                             (25,-20),(24,-25)])
    if data=='wrc':
        rain, hrtime, rlat, rlon, dates = openWRCstations()
        yrs = np.arange(1979,1999)
    elif data=='trmm':
        rain, hrtime, rlat, rlon, dates = openTRMMstations()
        yrs = np.arange(1998,2011)
    xy = np.hstack((rlon,rlat))
    ix = points_inside_poly(xy,summerregion)
    ix = np.where(ix)[0]
    rlat=rlat[ix];rlon=rlon[ix]
    #ixt = ixdtimes(dates,False,[10,11,12,01,02,03], False, False,mode='match')
    srain = sps.nanmean(rain[ix,:],axis=0)
    print srain.shape

    seasontot=[]
    for y in yrs:
        #md=[7,7, 1,31]
        ixt = ixdtimes(dates, [y,y+1], [10,03], [1,30],[0,0],mode='range')
        #print dates[ixt]
        #stotrain = my.sps.nansum(rain[ix,ixt])
        #print rain.shape,ix.shape,ixt.shape
        stotrain = np.nansum(srain[ixt])
        seasontot.append(stotrain)
    seasontot=np.asarray(seasontot)

    return seasontot, yrs

def monthlystationtotal(rain,dates,zscores=False,returnparams=False):
    '''Assumes rainfall data is station x time'''
    yrs = np.unique(dates[:-1,0])
    monthlytotal = np.zeros((rain.shape[0],len(yrs)*12))
    months = np.zeros((len(yrs)*12,4))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        for imn in xrange(12):
            mn = imn+1
            ix = np.where((dates[:,0]==yr) & (dates[:,1]==mn))[0]
            monthlytotal[:,iyr*12+imn] = np.sum(rain[:,ix],axis=1)
            #monthlytotal[:,iyr*12+imn] = rain[:,ix].sum(1)
            months[iyr*12+imn,:] = np.asarray([yr,mn,0,0])

    # Now calculate the anomalies
    meantotal = np.zeros((monthlytotal.shape[0],12))
    stdtotal = np.zeros((monthlytotal.shape[0],12))
    for i in xrange(12):
        ix = np.where(months[:,1] == i+1)[0]
        meantotal[:,i] = np.mean(monthlytotal[:,ix],axis=1)
        stdtotal[:,i]   = np.std(monthlytotal[:,ix],axis=1)

    anomtotal = monthlytotal - np.tile(meantotal,(1,len(yrs)))
    if zscores:
        anomtotal = (monthlytotal - np.tile(meantotal,(1,len(yrs))))/\
                                    np.tile(stdtotal,(1,len(yrs)))

    out=(monthlytotal, np.int32(months), anomtotal, meantotal)
    if returnparams:
        out=(monthlytotal, np.int32(months), anomtotal,meantotal,\
             np.tile(meantotal,(1,len(yrs))),np.tile(stdtotal,(1,len(yrs))),)
    return out

def seasonstationtotal(rain,dates,zscores=False,returnparams=False,\
                       season='coreseason'):
    '''Assumes rainfall data is station x time
       USAGE: season={'coreseason'} can also be
                      (str) 'fullseason'
                      (list of int) [11,12,1,2]
                      (tuple) (startmonth,startdate,seasonlen)
              This functionality actually just turns everything into a
              starting date plus season length which if not specified
              explicitly will be nmonths*30 days
    '''
    ### DEAL WITH DATES BEING [YR,MN,DY,HR] OR DATETIME OR DATETIME64
    if dates.ndim==2: nds=dates.shape[1]
    else: nds=0
    if (nds==4):
        datelist,yrslist=[],[]
        for d in dates:
            yr,mn,dy,hr = d
            datelist.append(np.datetime64(datetime.datetime(yr,mn,dy,hr)))
            yrslist.append(yr)
        dates=np.asarray(datelist)
    elif dates.dtype==np.dtype('datetime64[us]'):
        yrslist=[]
        for d in dates:
            yrslist.append(d.astype(np.ndarray).year)
    elif dates.dtype==np.dtype('O'):
        yrslist=[]
        for d in dates:
            yrslist.append(d.year)
    ### GET SEASON START AND LENGTH FOR MONTHS OR KEYWORD "SEASON"
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        seasonlength = "%d months" %(len(mns))
        startmonth, startdate = mns[0], 1
    elif isinstance(season,list):
        mns=season
        seasonlength = "%d months" %(len(mns))
        startmonth, startdate = mns[0], 1
    elif isinstance(season,tuple):
        startmonth,startdate,seasonlength = season
    seaslen, lenunit = seasonlength.split()
    if lenunit=='days': tdelta = datetime.timedelta(days=float(seaslen))
    elif lenunit=='weeks': tdelta = datetime.timedelta(weeks=float(seaslen))
    elif lenunit=='months': tdelta = datetime.timedelta(days=float(seaslen)*30)
    ### COMPUTE SEASONTOTALS
    yrs=np.unique(np.asarray(yrslist))
    seasontotal = np.zeros((rain.shape[0],len(yrs)))
    seasons = np.zeros((len(yrs)))
    for iyr, yr in enumerate(yrs):
        day1=datetime.datetime(yr,startmonth,startdate)
        daylast= day1 + tdelta
        ixt = ixdatetimes(dates,[day1,daylast],mode='range')
        seasons[iyr] = yr
        if isinstance(ixt,bool):
            print "Insufficient data for season starting:",yr
            seasontotal[:,iyr] = np.nan
            #seasons[iyr] = np.nan
            continue
        seasontotal[:,iyr] = np.sum(rain[:,ixt],axis=1)
    #irealt = np.where(~np.isnan(seasons))[0]
    #seasons = seasons[irealt]
    #seasontotal = seasontotal[:,irealt]
    ### NOW CALCULATE THE ANOMALIES
    meantotal = np.nanmean(seasontotal,axis=1)
    stdtotal        = np.nanstd(seasontotal,axis=1)
    anomtotal = seasontotal - np.tile(meantotal[:,np.newaxis],(1,len(seasons)))
    if zscores:
        anomtotal = (seasontotal - \
                     np.tile(meantotal[:,np.newaxis],(1,len(seasons))))\
                     /np.tile(stdtotal[:,np.newaxis],(1,len(seasons)))

    out=(seasontotal, np.int32(seasons), anomtotal, meantotal)
    if returnparams:
        out=(seasontotal, np.int32(seasons), anomtotal, meantotal,\
             np.tile(meantotal[:,np.newaxis],(1,len(yrs))),\
             np.tile(stdtotal[:,np.newaxis],(1,len(seasons))),)
    return out

def summerstations(raindata):
    '''Simply returns only the summer rainfall stations in a polygon defined
       within this method'''
    #summerregion=np.asarray([(24,-25),(28,-34),(35,-33),(33,-20),(25,-20),\
    #                         (24,-25)])
    summerregionsmaller=np.asarray([(23,-25),(28,-34),(35,-33),(35,-20),\
                             (25,-20),(23,-25)])
    summerregionbigger=np.asarray([(20,-25),(28,-34),(35,-33),(35,-15),\
                             (20,-15),(20,-25)])
    rain,dates,xy = raindata
    ix = points_inside_poly(xy,summerregionbigger)
    ix = np.where(ix)[0]
    summerrain = rain[ix,:]
    summerxy = xy[ix,:]

    return (summerrain,dates,summerxy), ix

def summerstations_regions(raindata,m=False):
    '''Simply returns only the summer rainfall stations as suggested in 
       Tennant & Hewitson 2001(? I think ?)

    USAGE: raindata (tuple) (rain,dates,xy)
    RETURNS: list of tuples (regionstations_data,dates,regionstations_lonlat)'''
    #summerregion=np.asarray([(24,-25),(28,-34),(35,-33),(33,-20),\
    #                         (25,-20),(24,-25)])
    summerregion=np.asarray([(23,-25),(28,-34),(35,-33),(35,-20),\
                             (25,-20),(23,-25)])
    summerregionsmaller=np.asarray([(23,-25),(28,-34),(35,-33),(35,-20),\
                             (25,-20),(23,-25)])
    summerregionbigger=np.asarray([(20,-25),(28,-34),(35,-33),(35,-15),\
                             (20,-15),(20,-25)])
    region1 = np.asarray([(29.,-22.5),(30.,-24.3),(32.1,-24.9),\
                          (32.,-22.),(29.,-22.5)])
    region2 = np.asarray([(26.5,-24.5),(29.,-27.),(29.,-29.7),(33.,-27.2),\
                          (32.1,-24.9),(30.,-24.3),(29.,-22.5),(26.5,-24.5)])
    region3 = np.asarray([(26.5,-24.5),(29.,-27.),(29.,-29.7),(29.,-30.),\
                          (26.5,-30.7),(23.5,-25.5),(26.5,-24.5)])
    region4 = np.asarray([(32.5,-27.5),(26.,-31.5),(28.,-34.),(33.,-28.),\
                          (32.5,-27.5)])
    if m:
        propstr='k--';lw=2.
        m.plot(region1[:,0],region1[:,1],propstr,lw=lw)
        plt.text(region1[:,0].mean(),region1[:,1].mean(),'1',\
                 fontsize=20,fontweight='bold')
        m.plot(region2[:,0],region2[:,1],propstr,lw=lw)
        plt.text(region2[:,0].mean(),region2[:,1].mean(),'2',\
                 fontsize=20,fontweight='bold')
        m.plot(region3[:,0],region3[:,1],propstr,lw=lw)
        plt.text(region3[:,0].mean(),region3[:,1].mean(),'3',\
                 fontsize=20,fontweight='bold')
        m.plot(region4[:,0],region4[:,1],propstr,lw=lw)
        plt.text(region4[:,0].mean(),region4[:,1].mean(),'4',\
                 fontsize=20,fontweight='bold')
        return


    rain,dates,xy = raindata
    regionrain=[]
    for summerregion in [region1, region2, region3, region4]:
        ix = points_inside_poly(xy,summerregion)
        ix = np.where(ix)[0]
        regrain = rain[ix,:]
        regxy = xy[ix,:]
        regionrain.append((regrain,dates,regxy))

    return regionrain

def TH01_regionrainindices(raindata,season=[11,12,1,2]):
    '''Computes normalized index for each of 4 regions described in 
       Tennant & Hewitson 2001

    USAGE: raindata (tuple) (rain,dates,xy)
    RETURNS: list of regional indices, timecoords'''
    regionrain=summerstations_regions(raindata)
    regionindices=[]
    for region in regionrain:
        rain, rdtime, lonlat = region
        totals, stime, anoms,climatol = \
                     seasonstationtotal(rain,rdtime,zscores=True,season=season)
        # Since compute normalized anomalies (zscores=True) 
        # just average to get average regional anomaly
        region_index = np.nanmean(anoms,axis=0)
        regionindices.append(region_index)
    return regionindices, stime

def readincutoffs(fname=\
                  '/home/neil/phd/cutofflowtracks_alicefavre/ncep2tracks.txt'):
    '''Reads in Alice Favre's cut-off low tracks: 
    A. Favre et al (2011) Clim. Dyn.  DOI:10.1007/s00382-011-1030-4'''
    fname='/home/neil/phd/cutofflowtracks_alicefavre/tracksNCEPDOE19792007'
    f = open(fname,'r')
    cnt=0
    tracks=[]
    for line in f:
        l=line.strip();l=l.strip('c(');l=l.strip(')')
        l = l.split(',')
        l[-1] = l[-1].strip()
        incr = len(l)/4
        track = []
        for i in xrange(incr):
            e = l[i::incr]
            e.append(cnt)
            entry= tuple(e)
            track.append(entry[::-1])
        tracks.extend(track)
        cnt+=1
    f.close()
    tracks = np.float32(tracks)
    # MAKE DATES
    start=[1979,1,1,0]
    doy=[];doyl=[]
    cnt=1
    for ed in monthends:
        tuppair = [(cnt,d) for d in range(1,ed+1)]
        doy.extend(tuppair);cnt+=1
    cnt=1
    for ed in monthends_leap:
        tuppair = [(cnt,d) for d in range(1,ed+1)]
        doyl.extend(tuppair);cnt+=1
    timekey=tracks[:,1]
    dates=[];yr=1978 # actually starts in 1979 but this year may need changing to a because of yr incrementing method (if timekey[i-1]-dky > 340: yr+=1)
    for i in xrange(len(timekey)):
        dky = int(timekey[i])
        if timekey[i-1]-dky > 340: yr+=1
        if yr%4 != 0: mn,dy=doy[dky-1]
        else: mn,dy=doyl[dky-1]
        dates.append((yr,mn,dy,0))

    dates=np.asarray(dates)

    return tracks, dates2hrnum(dates), dates

def readin_rwb(wb='AWB',pv=1.5,isk=340,\
               fdir='/home/neil/sogehome/data/rwb/fromThando/'):
    '''Read in Rossby Wave Breaking events from Thando Ndarana's 
       PV turn-over identifications
       Ndarana & Waugh (2011) J. Atmos. Sci.
    entries from track are: date,date.event_tag,ngridptsinoverturnregion,lon,
                            lat,PVcontour,AWB or CWB (1 or-1),isk'''
    class rwbtrack():
        def __init__(self,trk):
            self.lon, self.lat = trk[3], trk[4]
            self.nptsoverturned = trk[2]
            self.event_tag = trk[1]
            self.date = trk[0]
            self.PV, self.isentrope = trk[5], trk[7]
            breakflag=trk[6]
            if breakflag==1:
                self.AWB, self.CWB = True, False
            if breakflag==-1:
                self.AWB, self.CWB = False, True

    fname=fdir+wb+'.'+str(pv)+'.'+str(isk)+'.DJF'
    f=open(fname,'r')
    tracks=[];dates=[]
    for line in f:
        l=line.strip()
        l=l.split()
        entry= tuple(l)
        tracks.append(entry)
        yr,m,d=entry[0][:4],entry[0][4:6],entry[0][6:]
        dates.append((yr,m,d,0))
    f.close()

    tracks,dates=np.float32(tracks),np.int32(dates)
    trklist=[]
    for trk in tracks:
        trkobj = rwbtrack(trk)
        trklist.append(trkobj)
    return trklist, dates2hrnum(dates), dates

def readin_year_by_month(fname,missingval=-999.90):
    '''Read in indices data that have year by month format'''
    f=open(fname,'r')
    vals,dtime=[],[]
    for line in f:
        l=line.strip()
        yr,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12 = l.split()
        dts=[(yr,mn,1,0) for mn in range(1,13)]
        dtime.extend(dts)
        vals.extend([m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12])
    f.close()

    vals,dtime=np.float32(vals),np.int32(dtime)

    vals=np.where(vals==missingval,np.nan,vals)

    return vals,dates2hrnum(dtime),dtime


def mjo_phases(fname='/home/neil/data/indices/mjo/RMM1RMM2_WH04.txt'):

    f=open(fname)

    dtime=[];mjoix=[]
    for line in f:
        l=line.strip()
        year, month, day, RMM1, RMM2, phase, amplitude, forget = l.split()
        dtime.append((year,month,day,'0'))

        mjoix.append((RMM1,RMM2,phase,amplitude))
    f.close()

    dtime=np.int32(dtime)
    mjoix=np.float32(mjoix)

    return mjoix, dates2hrnum(dtime), dtime

def qbo_index(fname=\
            '/home/neil/sogehome/data/climate_indices/noaacpc.qbo.u50.index',\
              whichseries='standardized'):
    '''qboindex, time, dtime = qbo_index(fname="name file")
    '''
    datasrc = fname.split('/')[-1].split('.')[0]
    f = open(fname)
    if datasrc=='noaacpc':
        years=[]
        zonalwind=[]
        anomalies=[]
        zonalwind_std=[]
        datatype_switch=0
        for line in f:
            l=line.strip()
            linelist=l.split()
            ### Check if elements of line indicate not within datastream
            if len(linelist)<13: continue
            ### Increase dataswitch each time line including text "YEAR" occurs
            if linelist[0]=='YEAR':
                datatype_switch+=1
            else:
                linelist=np.float32(linelist)
            ### When first line element is a year (1979, etc.) add to data lists
            if isinstance(linelist,np.ndarray) and linelist[0]>1900.:
                if datatype_switch==1:
                    zonalwind.extend(linelist[1:])
                    years.append(linelist[0]) ### Same for other data too
                elif datatype_switch==2:
                    anomalies.extend(linelist[1:])
                elif datatype_switch==3:
                    zonalwind_std.extend(linelist[1:])
        f.close()  
        dtime=[]
        for yr in years:
            for mn in range(1,13):
                dtime.append(datetime.datetime(yr,mn,15,0))
        dtime = np.asarray(dtime)
        if whichseries=='standardized':
            qboindex = np.asarray(zonalwind_std)
        elif whichseries=='anomalies':
            qboindex = np.asarray(anomalies)
        elif whichseries=='original':
            qboindex = np.asarray(zonalwind)
        time = date2num(dtime,units="hours since 1800-01-01 00:00")
   
    return qboindex, time, dtime

def sst_indices(whichone='NINO3.4',fname=\
            '/home/neil/sogehome/data/climate_indices/noaacpc.sstoi.indices',\
                whichseries='standardized'):
    
    datasrc = fname.split('/')[-1].split('.')[0]
    f = open(fname,'r')
    if datasrc=='noaacpc':
        header=f.readline()
        f.close()
        columns=header.split()
        data=np.loadtxt(fname,skiprows=1)
        ichoose=columns.index(whichone)
        datachosen = data[:,ichoose]
        years, months = np.int32(data[:,0]), np.int32(data[:,1])
        dtime=[]
        for i,yr in enumerate(years):
            dtime.append(datetime.datetime(yr,months[i],15,0))
        dtime = np.asarray(dtime)
        if whichseries=='standardized':
            index=np.zeros((len(datachosen),))
            for i in xrange(12):
                ix = np.where(months[:] == i+1)
                index[ix] = sps.zscore(datachosen[ix])
        elif whichseries=='anomalies':
            index=np.zeros((len(datachosen),))
            for i in xrange(12):
                ix = np.where(months[:] == i+1)
                index[ix] = datachosen[ix]-datachosen[ix].mean()
        elif whichseries=='original':
            index = datachosen
        time = date2num(dtime,units="hours since 1800-01-01 00:00")
    
    return index, time, dtime

        

def monthlygrdvarmean(var,dates,zscores=False,returnparams=False):
    '''Assumes variable is gridded data is time x lat x lon'''
    yrs = np.unique(dates[:-1,0])
    monthlymean = np.zeros((len(yrs)*12,var.shape[1],var.shape[2]))
    months = np.zeros((len(yrs)*12,4))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        for imn in xrange(12):
            mn = imn+1
            ix = np.where((dates[:,0]==yr) & (dates[:,1]==mn))[0]
            if len(ix)==0   :monthlymean[iyr*12+imn,:] = np.nan
            else:monthlymean[iyr*12+imn,:] = sps.nanmean(var[ix,:],axis=0)
            months[iyr*12+imn,:] = np.asarray([yr,mn,0,0])

    # Now calculate the anomalies
    meantotal = np.zeros((12,monthlymean.shape[1],monthlymean.shape[2]))
    stdtotal = np.zeros((12,monthlymean.shape[1],monthlymean.shape[2]))
    for i in xrange(12):
        ix = np.where(months[:,1] == i+1)[0]
        meantotal[i,:] = sps.nanmean(monthlymean[ix,:],axis=0)
        stdtotal[i,:]   = sps.nanmean(monthlymean[ix,:],axis=0)

    anomtotal = monthlymean - np.tile(meantotal,(len(yrs),1,1))
    if zscores:
        anomtotal = (monthlymean - np.tile(meantotal,(len(yrs),1,1)))/\
                                   np.tile(stdtotal,(len(yrs),1,1))
    out=(monthlymean, np.int32(months), anomtotal, meantotal)
    if returnparams:
        out=(monthlymean, np.int32(months), anomtotal, meantotal,\
            np.tile(meantotal,(len(yrs),1,1)),np.tile(stdtotal,(len(yrs),1,1)),)
    return out

def seasongrdvarmean(var,dates,notconsecutiveyears=False,zscores=False,\
                     returnparams=False,season=[11,12,1,2]):
    '''Assumes variable is gridded data is time x lat x lon
       and that time is an np.array nt x 4 (year, mn, dy, hr)
       if however an array of datetime objects is given this is 
       converted to achieve this time construct
    '''
    if isinstance(dates[0],datetime.datetime):
        dates = np.asarray([(t.year, t.month, t.day, t.hour) for t in dates])
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season
    mst,mend=mns[0],mns[-1]
    dst,dend=1,monthends[mend]
    yrs = np.unique(dates[:-1,0])[:-1]
    if isinstance(notconsecutiveyears,np.ndarray):yrs=notconsecutiveyears
    seasonmean = np.zeros((len(yrs),var.shape[1],var.shape[2]))
    seasons = np.zeros((len(yrs)))
    for iyr in xrange(len(yrs)):
        yr=yrs[iyr]
        ixt = ixdtimes(dates,[yr,yr+1],[mst,mend],[dst,dend],[0,0],mode='range')
        seasonmean[iyr,:] = sps.nanmean(var[ixt,:],axis=0)
        seasons[iyr] = yr
    # Now calculate the anomalies
    meantotal = sps.nanmean(seasonmean,axis=0)
    stdtotal  = sps.nanstd(seasonmean,axis=0)

    anomtotal = seasonmean - np.tile(meantotal,(len(seasons),1,1))
    if zscores: anomtotal=(seasonmean -np.tile(meantotal,(len(seasons),1,1)))\
                                      /np.tile(meanstd,(len(seasons),1,1))

    out=(seasonmean, np.int32(seasons), anomtotal, meantotal)
    if returnparams:
        out=(seasonmean, np.int32(seasons), anomtotal,\
             meantotal,np.tile(meantotal,(len(seasons),1,1)),\
             np.tile(stdtotal,(len(seasons),1,1)),)
    return out

def seasongrdvarmean_deltat(var,dates,seasonstartdate,seasonlength="120 days",\
             notconsecutiveyears=False,zscores=False,returnparams=False):
    '''seasonmeans, seasons, anomalies, seasonalmean =\
       seasongrdvarmean_deltat(var,dates,seasonstartdate,\
                               seasonlength="120 days",\
                               notconsecutiveyears=False,zscores=False,\
                               returnparams=False)

      USAGE: var is your data nt x lat x lon
             dates is an array of datetime.datetime objects
             seasonstartdate is a tuple (month, day)
    '''
    startmonth,startdate = seasonstartdate
    seaslen, lenunit = seasonlength.split()
    if lenunit=='days': tdelta = datetime.timedelta(days=float(seaslen))
    elif lenunit=='weeks': tdelta = datetime.timedelta(weeks=float(seaslen))
    elif lenunit=='months':
        print 'Converting to days, 30 days per month'
        tdelta = datetime.timedelta(days=float(seaslen)*30)
    yrslist=np.asarray([t.year for t in dates])
    yrs=np.unique(yrslist)
    if isinstance(notconsecutiveyears,np.ndarray):yrs=notconsecutiveyears
    seasonmean = np.zeros((len(yrs),var.shape[1],var.shape[2]))
    seasons = np.zeros((len(yrs)))
    nanfield = np.zeros((var.shape[1],var.shape[2]))
    nanfield[:] = np.nan
    for iyr, yr in enumerate(yrs):
        day1=datetime.datetime(yr,startmonth,startdate)
        daylast= day1 + tdelta
        ixt = ixdatetimes(dates,[day1,daylast],mode='range')
        seasons[iyr] = yr
        if isinstance(ixt,bool):
            print "Insufficient data for season starting:",yr
            seasonmean[iyr,...] = nanfield
            #seasons[iyr] = np.nan
            continue
        seasonmean[iyr,:] = sps.nanmean(var[ixt,:],axis=0)
    #irealt = np.where(~np.isnan(seasons))[0]
    #seasons = seasons[irealt]
    #seasonmean = seasonmean[irealt,...]
     
    # Now calculate the anomalies
    meantotal = sps.nanmean(seasonmean,axis=0)
    stdtotal  = sps.nanstd(seasonmean,axis=0)

    anomtotal = seasonmean - np.tile(meantotal,(len(seasons),1,1))
    if zscores: anomtotal=(seasonmean -np.tile(meantotal,(len(seasons),1,1)))\
                                      /np.tile(meanstd,(len(seasons),1,1))
    out=(seasonmean, np.int32(seasons), anomtotal, meantotal)
    if returnparams:
        out=(seasonmean, np.int32(seasons), anomtotal,\
             meantotal,np.tile(meantotal,(len(seasons),1,1)),\
             np.tile(stdtotal,(len(seasons),1,1)),)
    return out

def seasongrdvarstd_deltat(var,dates,seasonstartdate,seasonlength="120 days",\
             notconsecutiveyears=False,zscores=False,returnparams=False):
    '''seasonstds, seasons, anomalies, seasonalstd =\
       seasongrdvarmean_deltat(var,dates,seasonstartdate,\
                               seasonlength="120 days",
                               notconsecutiveyears=False,zscores=False,\
                               returnparams=False)

      USAGE: var is your data nt x lat x lon
             dates is an array of datetime.datetime objects
             seasonstartdate is a tuple (month, day)
    '''
    startmonth,startdate = seasonstartdate
    seaslen, lenunit = seasonlength.split()
    if lenunit=='days': tdelta = datetime.timedelta(days=float(seaslen))
    elif lenunit=='weeks': tdelta = datetime.timedelta(weeks=float(seaslen))
    elif lenunit=='months':
        print 'Converting to days, 30 days per month'
        tdelta = datetime.timedelta(days=float(seaslen)*30)
    yrslist=np.asarray([t.year for t in dates])
    yrs=np.unique(yrslist)
    if isinstance(notconsecutiveyears,np.ndarray):yrs=notconsecutiveyears
    seasonmean = np.zeros((len(yrs),var.shape[1],var.shape[2]))
    seasons = np.zeros((len(yrs)))
    nanfield = np.zeros((var.shape[1],var.shape[2]))
    nanfield[:] = np.nan
    for iyr, yr in enumerate(yrs):
        day1=datetime.datetime(yr,startmonth,startdate)
        daylast= day1 + tdelta
        ixt = ixdatetimes(dates,[day1,daylast],mode='range')
        seasons[iyr] = yr
        if isinstance(ixt,bool):
            print "Insufficient data for season starting:",yr
            seasonmean[iyr,...] = nanfield
            #seasons[iyr] = np.nan
            continue
        seasonmean[iyr,:] = sps.nanstd(var[ixt,:],axis=0)
    #irealt = np.where(~np.isnan(seasons))[0]
    #seasons = seasons[irealt]
    #seasonmean = seasonmean[irealt,...]
     
    # Now calculate the anomalies
    meantotal = sps.nanmean(seasonmean,axis=0)
    stdtotal  = sps.nanstd(seasonmean,axis=0)

    anomtotal = seasonmean - np.tile(meantotal,(len(seasons),1,1))
    if zscores: anomtotal=(seasonmean -np.tile(meantotal,(len(seasons),1,1)))\
                                      /np.tile(meanstd,(len(seasons),1,1))
    out=(seasonmean, np.int32(seasons), anomtotal, meantotal)
    if returnparams:
        out=(seasonmean, np.int32(seasons), anomtotal,\
             meantotal,np.tile(meantotal,(len(seasons),1,1)),\
             np.tile(stdtotal,(len(seasons),1,1)),)
    return out

def seasoncomposite(data,years,subsetindices,ttest=True,confidence=.95,\
    twosample=False):
    ''' compositemean, compositeyears, yrs2comp_correlation = 
              seasoncomposite(data,years,subsetindices,ttest=True)
        
        USAGE: data (array) with time as first dimension
                            NOTE: data can be tuple(u,v) for vector fields
                                  and a Hoteling test will be used instead
                                  but there will be no spatial correlation
                                  returned
               years (array) specifying the season year of each timestep in data
               subsetindices are use compute the composite mean
               ttest={True} however this 1-sample t-test will only be correct 
                            if data is normally distributed.
               confidence={0.95} is the confidence interval for t-test
        
        Returns: compositemean which will be a masked array if ttest=True
                 compositeyears
                 yrs2comp_correlation which gives the spatial correlation of
                                      each compositeyear with the compositemean
    '''
    ixsub = subsetindices
    if isinstance(data,np.ndarray):
        meancomp = np.nanmean(data[ixsub,...],axis=0)
        if ttest and twosample:
            meancomp,pval=ttestmask_2samp(data,data[ixsub,...],\
                          confidence=confidence)
        elif ttest and not twosample:
            meancomp,pval=ttestmask(data,meancomp,confidence=confidence)
        r = np.zeros(ixsub.shape)
        for i, ix in enumerate(ixsub):
            corr = np.corrcoef(meancomp.ravel(),data[ix,...].ravel())
            r[i] = corr[0,1]
    elif isinstance(data,tuple):
        ### THIS IS NOT TESTED AND SEEMS TO NOT WORK PROPERLY!
        u, v = data
        meancompu = np.nanmean(u[ixsub,...],axis=0)
        meancompv = np.nanmean(v[ixsub,...],axis=0)
        if ttest:
            meancompu,meancompv,pval=hoteltestmask_1samp((u,v),\
                           (meancompu,meancompv),confidence=confidence)
            meancomp = (meancompu, meancompv)
        ### Not sure how best to compute this correlation since it is u and v
        r=False
        #r = np.zeros(ixsub.shape)
        #for i, ix in enumerate(ixsub):
        #    corr = np.corrcoef(meancomp.ravel(),data[ix,...].ravel())
        #    r[i] = corr[0,1]

    return meancomp, years[subsetindices], r
        

def ttestmask(dataobs,datacomp,confidence=.95):
    '''Test data sample against another data  mean (eg. test confidence of 
       composite not belonging to climatology distribution)'''
    ### Remove all nan from dataobs
    ndm = dataobs.ndim
    if ndm==3:
        irealt = np.where(~np.isnan(dataobs[:,0,0]))[0]
        dataobsnonan = dataobs[irealt,...]
    elif ndm==2:
        irealt = np.where(~np.isnan(dataobs[:,0]))[0]
        dataobsnonan = dataobs[irealt,...]
    if ndm==1:
        irealt = np.where(~np.isnan(dataobs))[0]
        dataobsnonan = dataobs[irealt]
    ### Now t-test composite against obs sample
    ttstat, pval = sps.ttest_1samp(dataobsnonan,datacomp,axis=0)
    pvalmsk = pval > (1-confidence)
    #pv=np.where(pvalmsk,np.nan,pval)
    return np.ma.MaskedArray(datacomp,mask=pvalmsk,dtype=np.float32), pval

def ttestmask_2samp(dataobs,datasamp,confidence=.95):
    '''Test data sample distribution against another data sample distribution
              (eg. test confidence of on sample not belonging within other
                   sample distribution)'''
    ### Remove all nan from dataobs and datasamp samples
    ndm = dataobs.ndim
    if ndm==3:
        irealt = np.where(~np.isnan(dataobs[:,0,0]))[0]
        dataobsnonan = dataobs[irealt,...]
    elif ndm==2:
        irealt = np.where(~np.isnan(dataobs[:,0]))[0]
        dataobsnonan = dataobs[irealt,...]
    if ndm==1:
        irealt = np.where(~np.isnan(dataobs))[0]
        dataobsnonan = dataobs[irealt]
    ndm = datasamp.ndim
    if ndm==3:
        irealt = np.where(~np.isnan(datasamp[:,0,0]))[0]
        datasampnonan = datasamp[irealt,...]
    elif ndm==2:
        irealt = np.where(~np.isnan(datasamp[:,0]))[0]
        datasampnonan = datasamp[irealt,...]
    if ndm==1:
        irealt = np.where(~np.isnan(datasamp))[0]
        datasampnonan = datasamp[irealt]
    ### Now t-test composite against obs sample
    ttstat, pval = sps.ttest_ind(dataobsnonan,datasampnonan,axis=0)
    pvalmsk = pval > (1-confidence)
    #pvalmsk=np.tile(pvalmsk,(datasamp.shape[0],1,1))
    #pv=np.where(pvalmsk,np.nan,pval)
    datacomp=datasamp.mean(0)
    return np.ma.MaskedArray(datacomp,mask=pvalmsk,dtype=np.float32), pval

def hoteltestmask(datareturn_uv,datacomp_uv,confidence=.95):
    '''Implemented as 2 independant sample version
       see hoteltestmask_1samp for 1 sample version

    USAGE: datareturn_uv is tuple(u,v) where u is time x lat x lon
           datareturn also denotes this will be return as masked array
           datacomp_uv has same dimensions but is the comparison sample 
                      (eg. climatology run)

    RETURNS: u, v, pvals

    NOTE: Built mainly from information at 
          http://en.wikipedia.org/wiki/Hotelling%27s_T-square_distribution
          based on Hotelling's original paper and 
          Mardia, Kent and Biddy (1979) text book
    '''
    ### Note that at present this test does not have good times with nans in the input
    nt,nx,ny=datareturn_uv[0].shape;
    nv=2.
    Xu=datareturn_uv[0].reshape(nt,nx*ny,order='F')
    Xv=datareturn_uv[1].reshape(nt,nx*ny,order='F')
    Yu=datacomp_uv[0].reshape(nt,nx*ny,order='F')
    Yv=datacomp_uv[1].reshape(nt,nx*ny,order='F')

    pvals=np.ones((1,nx*ny))
    for i in xrange(nx*ny):
        X=np.vstack((Xu[:,i],Xv[:,i])).T
        Y=np.vstack((Yu[:,i],Yv[:,i])).T
        Xmn,Ymn=X.mean(0),Y.mean(0)
        dXYmn=np.asmatrix(Xmn-Ymn)
        Xan,Yan=X-np.tile(Xmn,(nt,1)), Y-np.tile(Ymn,(nt,1))

        W = np.linalg.inv(np.asmatrix(np.cov((Xan-Yan).T)))
        T2=(nt*nt/(nt+nt))*dXYmn*W*dXYmn.H
        F=(2.*nt-nv-1.)/(2.*nt-2.)/nv*T2
        pvals[0,i] = sps.f.cdf(F, nv, nt-nv)

    pvalmsk=pvals > (1-confidence)
    pvalmsk = np.tile(pvalmsk,(nt,1))
    pvalmsk = pvalmsk.reshape(nt,nx,ny,order='F')

    u = np.ma.MaskedArray(datareturn_uv[0],mask=pvalmsk,dtype=np.float32)
    v = np.ma.MaskedArray(datareturn_uv[1],mask=pvalmsk,dtype=np.float32)

    u,v = u.mean(0),v.mean(0)

    return u,v, pvals.reshape(nx,ny,order='F')

def hoteltestmask_1samp(dataobs_uv,datacomp_uv,confidence=.95):
    '''Test U,V sample against another U,V mean (eg. test confidence of 
       composite U,V not belonging to climatology distribution)

    USAGE: dataobs_uv is tuple(u,v) where u is time x lat x lon
           datacomp_uv has same dimensions but is the comparison mean composite 

    RETURNS: u, v, pvals

    NOTE: Built mainly from information at 
          http://en.wikipedia.org/wiki/Hotelling%27s_T-square_distribution
          based on Hotelling's original paper and 
          Mardia, Kent and Biddy (1979) text book'''
    ### Note that at present this test does not have good times with nans in the input
    nt,nx,ny=dataobs_uv[0].shape
    nv=2.
    mu=datacomp_uv[0].reshape(nx*ny,order='F')
    mv=datacomp_uv[1].reshape(nx*ny,order='F')
    Xu=dataobs_uv[0].reshape(nt,nx*ny,order='F')
    Xv=dataobs_uv[1].reshape(nt,nx*ny,order='F')

    pvals=np.ones((1,nx*ny))
    for i in xrange(nx*ny):
        m=np.vstack((mu[i],mv[i])).T
        X=np.vstack((Xu[:,i],Xv[:,i])).T
        Xmn=X.mean(0)
        dXmn=np.asmatrix(Xmn-m)
        Xan=X-np.tile(Xmn[np.newaxis,:],(nt,1))

        W = np.linalg.inv(np.asmatrix(np.cov(Xan.T)))
        T2=nt*dXmn*W*dXmn.H
        F=(nt-nv)/nv/(nt-1.)/nv*T2
        pvals[0,i] = sps.f.cdf(F, nv, nt-nv)

    pvalmsk=pvals > (1-confidence)
    pvalmsk = pvalmsk.reshape(nx,ny,order='F')

    u = np.ma.MaskedArray(datacomp_uv[0],mask=pvalmsk,dtype=np.float32)
    v = np.ma.MaskedArray(datacomp_uv[1],mask=pvalmsk,dtype=np.float32)

    return u,v, pvals

def true_gcoords(longitudes,latitudes,latp=37.5,lonp=177.5,meshgrid=False):
    '''geolat,geolon = mytools.truegcoords(longitudes,latitudes,\
                                     latp=37.5,lonp=177.5,meshgrid=False)
    
    Returns the true latitude, longitude values:

    If meshgrid=True this returns meshgrid as expected 
    '''
    lambdap = np.deg2rad(lonp)-np.pi
    phip = np.pi/2. - np.deg2rad(latp)
    if meshgrid: 
        xx,yy = np.meshgrid(longitudes,latitudes)
        rlambda,rphi =xx.flatten(),yy.flatten()
    else: 
        rlambda,rphi = longitudes,latitudes
    rlambda = np.deg2rad(rlambda)
    rphi = np.deg2rad(rphi);
    tlam = (np.cos(phip)*np.sin(lambdap)*np.cos(rphi)*np.cos(rlambda) + \
    np.cos(lambdap)*np.cos(rphi)*np.sin(rlambda) - \
    np.sin(phip)*np.sin(lambdap)*np.sin(rphi)) / \
    (np.cos(phip)*np.cos(lambdap)*np.cos(rphi)*np.cos(rlambda) - \
    np.sin(lambdap)*np.cos(rphi)*np.sin(rlambda) - \
    np.sin(phip)*np.cos(lambdap)*np.sin(rphi))
    lmbda = np.rad2deg(np.arctan(tlam))
    sphi = np.sin(phip)*np.cos(rphi)*np.cos(rlambda) + np.cos(phip)*np.sin(rphi)
    phi = np.rad2deg(np.arcsin(sphi))
    if meshgrid: lmbda,phi= lmbda.reshape(xx.shape),phi.reshape(xx.shape)
    return lmbda, phi

def plotshadedline(x,y,z,alphafactor=False,zlim=False,colormap='jet',\
                   alpha=0.5,lw=1):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=plt.get_cmap(colormap))
    if zlim:
        lc = LineCollection(segments, cmap=plt.get_cmap(colormap),\
        norm=plt.Normalize(zlim[0],zlim[1],clip=True),alpha=alpha)
    lc.set_array(z)
    lc.set_linewidth(lw)
    plt.gca().add_collection(lc)
    if alphafactor:
        plt.show()
        ecols=lc.get_edgecolors()
        ecols[:,3]=alphafactor
        lc.set_edgecolors(ecols)
        plt.cla()
        plt.gca().add_collection(lc)
    return lc

def xunit_as_pts(xunit,f=False,ax=False):
    if not f and not ax:
        f=plt.gcf();ax=plt.gca()
    dpi=f.dpi
    fwidth_inches = f.get_figwidth()
    figdotswidth = fwidth_inches*dpi
    ### Axes width
    axbox = ax.get_position()
    awidthfraction = bbox.width
    axdotswidth = figdotswidth*awidthfraction
    xunits = ax.get_xlim()
    xwidth = np.abs(xunits[1]-xunits[0])
    nxunits = xwidth/xunit
    pts_of_xunit = axdotswidth/nxunits
    return pts_of_xunit

def yunit_as_pts(yunit,f=False,ax=False):
    if not f and not ax:
        f=plt.gcf();ax=plt.gca()
    dpi=f.dpi
    fheight_inches = f.get_figheight()
    figdotsheight = fheight_inches*dpi
    ### Axes height
    axbox = ax.get_position()
    aheightfraction = bbox.height
    axdotsheight = figdotsheight*aheightfraction
    yunits = ax.get_ylim()
    yheight = np.abs(yunits[1]-yunits[0])
    nyunits = yheight/yunit
    pts_of_yunit = axdotsheight/nyunits
    return pts_of_yunit

def plotmaprectangle(m,rect,rectype='corners',geocoords=True):
    if rectype=='corners':
        ln1,lt1=rect[0];ln2,lt2=rect[1];ln3,lt3=rect[2];ln4,lt4=rect[3];
    s1=np.asarray((np.linspace(ln1,ln2),np.linspace(lt1,lt2)))
    s2=np.asarray((np.linspace(ln2,ln3),np.linspace(lt2,lt3)))
    s3=np.asarray((np.linspace(ln3,ln4),np.linspace(lt3,lt4)))
    s4=np.asarray((np.linspace(ln4,ln1),np.linspace(lt4,lt1)))
    for i in [1,2,3,4]:
        exec('s%d[0,:],s%d[1,:]=true_gcoords(s%d[0,:],s%d[1,:])'%(i,i,i,i))
    sq=np.hstack((s1,s2,s3,s4))
    ### Project region on to map
    xs,ys=m(sq[0,:],sq[1,:])
    ### Deal with spurios values where line outside of map vision
    nanny=np.tile(np.nan,200)
    xs=np.where(xs>1.0e+15,nanny,xs )
    m.plot(xs,ys,lw=2.,color='k',ls=':')

def labelaxis():
    def labax(x, pos):
        'The two args are the value and tick position'
        return '$%1.1fM' % (x)
    return FuncFormatter(labax)

def closestlonlat(lon,lat,pt):
    '''Very simple function that return the lon and lat values and indices of
       that are closest to a given point, pt.
       USEFUL for find the grid point closest to a station'''
    ln,lt=pt
    londiff,latdiff = np.abs(lon-ln), np.abs(lat-lt)
    lndmin, ltdmin = londiff.min(), latdiff.min()
    iln = (londiff==lndmin).nonzero()[0]
    ilt = (latdiff==ltdmin).nonzero()[0]
    return lon, lat, iln, ilt

def spheredistance(lon1, lat1, lon2, lat2):
    '''Calculated great-circle distances in metres'''
    R=6367442.76
    # Determine proper longitudinal shift
    dlon=np.abs(lon2-lon1)
    dlon=np.where(dlon>=180,360-dlon,dlon)
    #  Convert Decimal degrees to radians.
    deg2rad=np.pi/180
    lat1=lat1*deg2rad
    lat2=lat2*deg2rad
    dlon=dlon*deg2rad
    #
    #  Compute the distances
    t1 = np.sin(dlon) * np.cos(lat2)
    t2 = np.sin(lat2) * np.cos(lat1)
    t3 = np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    dist=R*np.arcsin(np.sqrt(t1**2 + (t2-t3)**2))

    return dist

def xlim(l1,l2,fontweight='bold',fontsize=13.):
    plt.xlim(l1,l2)
    ax=plt.gca()
    labels=[l for l in ax.xaxis.get_ticklabels()]
    for lb in labels:
        lb.set_fontsize(fontsize)
        lb.set_fontweight(fontweight)

def ylim(l1,l2,fontweight='bold',fontsize=13.):
    plt.ylim(l1,l2)
    ax=plt.gca()
    labels=[l for l in ax.yaxis.get_ticklabels()]
    for lb in labels:
        lb.set_fontsize(fontsize)
        lb.set_fontweight(fontweight)
def ytickfonts(fontweight='bold',fontsize=13.,color='k',rotation=False):
    ax=plt.gca()
    labels=[l for l in ax.yaxis.get_ticklabels()]
    for lb in labels:
        lb.set_fontsize(fontsize)
        lb.set_fontweight(fontweight)
        lb.set_color(color)
        if rotation: lb.set_rotation(rotation)

def axistick2nd(axis,tickarr,scaling=1.):
    labels=[str(d) for d in tickarr[1::2]*scaling]
    labels=','+(',,').join(labels)+','
    labels=labels.split(',')
    if len(labels)>len(tickarr):labels=labels[:-1]
    exec('plt.'+axis+'ticks(tickarr,labels)')

def xtickfonts(fontweight='bold',fontsize=13.,color='k',rotation=False):
    ax=plt.gca()
    labels=[l for l in ax.xaxis.get_ticklabels()]
    for lb in labels:
        lb.set_fontsize(fontsize)
        lb.set_fontweight(fontweight)
        lb.set_color(color)
        if rotation: lb.set_rotation(rotation)

# Basic function for mkdir -p functionality
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise