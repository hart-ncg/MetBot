# MetBlobs wrapper
#
# Format to follow for all variables
# DSET-VAR-LEV-DERIVE
# DSET-VAR-LEV-DERIVE-{EXP}{ENS} (for flavours experiments [maybe ccam ouput too])
# When time subsets of dsets are used, this should be denoted
# DSET-VAR-LEV-DERIVE-{subsetdescription}

import numpy as np
import mynetcdf as mync
import mytools as my
import MetBlobs as blb
import SynopticAnatomy as sy
import SynopticPlot as syp
import datetime
import matplotlib.pyplot as plt
import cPickle
import time as tmr
import gc
tstart=tmr.time()
#import glob
#from mpl_toolkits.basemap import num2date

def datastretched(vrb,lat,lon,v,sub='SA'):
    '''This function copies first few steps of Metblobs()
    by taking data and then stretching it according to filterdict.py
    However it then uses this and the threshold value to calculate grid point frequencies for
    values that meet the threshold criteria'''
    dct=blb.blobfilters[sub+'cloudband'][v]
    domain=dct['ROI']
    dstretch=dct['stretch']
    dthresh, highlow = dct['thresh']
    data, thrs=blb.Stretch(vrb,dthresh,tval=dstretch)
    print "Stretched variable:",v
    if highlow=='low':
        #thrs=255-thrs
        datamsk = data<=thrs
    elif highlow=='high':
        datamsk = data>=thrs

    return ~datamsk

def rawdatamasks(varlist):

    noaaolr=False
    ncep2=False
    time6,time24=[],[]

    newvarlist=varlist[:]
    retlist=[]
    for v in varlist:
        dset,varstr, lev, drv = v.split('-')
        retlist.append('datam'+varstr+'_'+lev)
        if dset=='noaa' and varstr=='olr':
            noaaolr = True
            newvarlist.remove(v)
        if dset=='ncep2': ncep2=True
        if ncep2:
            if varstr=='uvmag':
                newvarlist.remove(v)
                newvarlist.extend(['ncep2-uwnd-'+lev+'-del2','ncep2-vwnd-'+lev+'-del2'])
            if varstr=='quvmag':
                newvarlist.remove(v)
                newvarlist.extend(['ncep2-qu-'+lev+'-del2','ncep2-qv-'+lev+'-del2'])


    sub="SA"
    picdir="/home/neil/work/computervision/metblobs/"+sub+"/pickles/"

    hgtthresh=1.3
    jetthresh=1.0
    moistthresh=0.5

    # CONSTRUCT A REFTIME FROM CENTRES OF THE SEASON 01-JAN
    refsubset=True
    hrwindow=24*30*3
    refdates=[]
    for yr in np.arange(1979,2012):
        refdates.append((yr,01,01,0))
    refdates=np.asarray(refdates)
    reftime = my.dates2hrnum(refdates)

    mycm = plt.cm.PuBu
    # OPENING FILES
    # Liebmann & Smith Interpolated OLR
    if noaaolr:
        v="noaa-olr-0-0"
        print '\n',v
        dset,varstr, lev, drv = v.split('-')
        dsrc="/home/neil/data/olr/"
        #olr, time, lat, lon, dtime = mync.openolr(dsrc+'olr.NDJF.Fauch08.nc','olr')
        olr, time, lat, lon, dtime = mync.openolr(dsrc+'olr.day.mean.nc','olr',subs=sub)
        time=time-15769752.0
        olr=olr[1675:,:,:];time=time[1675:];dtime=dtime[1675:] # gets all the data since 1979-01-01 00:00
        if refsubset: exec("ixt, [time, "+varstr+", dtime] = my.ixtwindow(reftime,time,hrwindow,time,"+varstr+",dtime)")
        time24=time
        exec('datamsk'+varstr+'_'+lev+' = datastretched('+varstr+',lat,lon,v,sub=sub)')
        exec('datam'+varstr+'_'+lev+' = np.ma.MaskedArray('+varstr+',mask=datamsk'+varstr+'_'+lev+',dtype=np.float32)')
        #m, f = blb.SAfrBasemap(lat[3:-6],lon[3:-3],drawstuff=True,prj='cyl',fno=1,rsltn='l')
        #exec('plt.figure();syp.redrawmap(m);m.pcolor(lon,lat,nfreq'+varstr+'_'+lev+',cmap=mycm);plt.clim(0.01,0.5);\
        #plt.title("'+varstr+'_'+lev+'")')

    # NCEP2 subset file of your choice (created in by ncep2subset in subsetfiles.py)
    if ncep2:
        print "Analysing ncep2 data set..."
        #hgtlist=['ncep2-hgt-850-del2']#,'ncep2-hgt-500-del2','ncep2-hgt-250-del2']
        #uvlist=['ncep2-uwnd-250-del2','ncep2-vwnd-250-del2']#,'ncep2-uwnd-500-del2','ncep2-vwnd-500-del2','ncep2-uwnd-250-del2','ncep2-vwnd-250-del2']
        #quvlist=['ncep2-qu-700-del2','ncep2-qv-700-del2','ncep2-qu-500-del2','ncep2-qv-500-del2','ncep2-qu-850-del2','ncep2-qv-850-del2']
        #quvlist=['ncep2-qu-850-del2','ncep2-qv-850-del2']
        #varlist=[]
        #varlist.extend(hgtlist)
        #varlist.extend(uvlist)
        #varlist.extend(quvlist)
        for v in newvarlist:
            print '\n',v
            # OPEN FILE
            sbst=sub
            dset, varstr, lev, deriv = v.split('-')
            dsrc="/home/neil/data/ncep2/6hourly/"+varstr+"/"
            ncfileout = dsrc+"%s.%s_%s.nc" %(varstr, sbst, lev)
            exec(varstr+", time, lat, lon, dtime = mync.open"+dset+"(ncfileout,\'"+varstr+"\')")
            # SUBSET IN TIME BY REFTIME
            if refsubset: exec("ixt, [time, "+varstr+", dtime] = my.ixtwindow(reftime,time,hrwindow,time,"+varstr+",dtime)")
            time6=time
            # CALCULATION OF DERIVED VARIABLES
            # Laplacian of pressure surface
            if deriv=='del2' and varstr=='hgt':
                print "Get laplacian of hgt..."
                del2 = my.del2(hgt,time,lat,lon)
                del2 = my.coriolosscale(del2,lat,lon)
                slims, thrs = blb.calclims(del2,scalestd=hgtthresh,disttype='normal')
                blb.blobfilters[sub+'cloudband'][v]['stretch'] = slims
                blb.blobfilters[sub+'cloudband'][v]['thresh'] = (thrs[1],'high')
                exec('datamsk'+varstr+'_'+lev+' = datastretched(del2,lat,lon,v,sub=sub)')
                exec('datam'+varstr+'_'+lev+' = np.ma.MaskedArray('+varstr+',mask=datamsk'+varstr+'_'+lev+',dtype=np.float32)')
                del del2, hgt
                #m, f = blb.SAfrBasemap(lat[3:-6],lon[3:-3],drawstuff=True,prj='cyl',fno=1,rsltn='l')
                #exec('plt.figure();syp.redrawmap(m);m.pcolor(lon,lat,nfreq'+varstr+'_'+lev+',cmap=mycm);plt.clim(0.01,0.9);\
                #plt.title("'+varstr+'_'+lev+'")')
                #plt.colorbar()

            # Laplacian of u,v magnitude
            elif 'uwnd' in locals() and 'vwnd' in locals() and deriv=='del2':
                print "Get laplacian of vector wind magnitude..."
                vmag=dset+"-uvmag-"+lev+"-"+deriv
                dset, varstr, lev, deriv = vmag.split('-')
                mag, degdir = my.magdir(uwnd,vwnd);del uwnd, vwnd
                del2 = my.del2(mag,time,lat,lon); del degdir
                del2 = my.coriolosscale(del2,lat,lon)
                slims, thrs = blb.calclims(del2,scalestd=jetthresh,disttype='normal')
                blb.blobfilters[sub+'cloudband'][vmag]['stretch'] = slims
                blb.blobfilters[sub+'cloudband'][vmag]['thresh'] = (thrs[0],'low')
                exec('datamsk'+varstr+'_'+lev+' = datastretched(del2,lat,lon,vmag,sub=sub)')
                exec('datam'+varstr+'_'+lev+' = np.ma.MaskedArray(mag,mask=datamsk'+varstr+'_'+lev+',dtype=np.float32)')
                #exec('datam'+varstr+'_'+lev+' = np.ma.MaskedArray(v,mask=datamsk'+varstr+'_'+lev+',dtype=np.float32)')
                del del2,mag
                #m, f = blb.SAfrBasemap(lat[3:-6],lon[3:-3],drawstuff=True,prj='cyl',fno=1,rsltn='l')
                #exec('plt.figure();syp.redrawmap(m);m.pcolor(lon,lat,nfreq'+varstr+'_'+lev+',cmap=mycm);plt.clim(0.01,0.5);\
                #plt.title("'+varstr+'_'+lev+'")')

            # Laplacian of qu,qv magnitude
            elif 'qu' in locals() and 'qv' in locals() and deriv=='del2':
                print "Get laplacian of vector moisture flux magnitude..."
                vmag=dset+"-quvmag-"+lev+"-"+deriv
                dset, varstr, lev, deriv = vmag.split('-')
                mag, degdir = my.magdir(qu,qv);del qu, qv
                del2 = my.del2(mag,time,lat,lon);del degdir
                #del2 = my.coriolosscale(del2,lat,lon)
                slims, thrs = blb.calclims(del2,scalestd=moistthresh,disttype='normal')
                blb.blobfilters[sub+'cloudband'][vmag]['stretch'] = slims
                blb.blobfilters[sub+'cloudband'][vmag]['thresh'] = (thrs[0],'low')
                exec('datamsk'+varstr+'_'+lev+' = datastretched(del2,lat,lon,vmag,sub=sub)')
                exec('datam'+varstr+'_'+lev+' = np.ma.MaskedArray(mag,mask=datamsk'+varstr+'_'+lev+',dtype=np.float32)')
                #exec('datam'+varstr+'_'+lev+' = np.ma.MaskedArray(v,mask=datamsk'+varstr+'_'+lev+',dtype=np.float32)')
                del del2,mag
                #m, f = blb.SAfrBasemap(lat[3:-6],lon[3:-3],drawstuff=True,prj='cyl',fno=1,rsltn='l')
                #exec('plt.figure();syp.redrawmap(m);m.pcolor(lon,lat,nfreq'+varstr+'_'+lev+',cmap=mycm);plt.clim(0.01,0.5);\
                #plt.title("'+varstr+'_'+lev+'")')

        exec('out=('+", ".join(retlist)+')')

        return out, time24, time6
