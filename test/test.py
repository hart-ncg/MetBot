# MetBlobs wrapper
#
# Format to follow for all variables
# DSET-VAR-LEV-DERIVE
# DSET-VAR-LEV-DERIVE-{EXP}{ENS} (for flavours experiments [maybe ccam ouput too])
# When time subsets of dsets are used, this should be denoted
# DSET-VAR-LEV-DERIVE-{subsetdescription}
import numpy as np
import datetime
import matplotlib.pyplot as plt
import cPickle
import time as tmr
import gc
import sys,os
cwd=os.getcwd()
### Add path to MetBot modules and import
sys.path.append(cwd+'/..')
import MetBot.mynetcdf as mync
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.SynopticAnatomy as sy
tstart=tmr.time()

### Ensure only look at Southern Africa
sub="SA"
### Place to store output files: CHANGE ACCORDINGLY
picdir=cwd+"/"
### Where to save the figs
figdir=picdir
noaaolr=True
noaaolrall=True  # This is for tracking pre and post TTT cloud blobs
noaaolrfull=True # This is for tracking pre and post TTT cloud blobs
testyear=True    # Only uses first 365 days of olr data
getdistr=True    # Save a figure showing histogram of OLR values
getmbs=True      # Actually run the MetBot algorithm
showblb=False    # Show the blobs while running
refsubset=True   # This is used if noaaolr=True to only look in time window
hrwindow=49      # ... close (49 hours/ 2days) to flagged cloud band days
synoptics=True   # Build tracks of cloud blobs that become TTT cloud bands
                 # ... which are then used to build TTT events.

# OPENING FILES
# Liebmann & Smith Interpolated OLR
if noaaolr:
    v="noaa-olr-0-0"
    dset,varstr, lev, drv = v.split('-')
    dsrc=picdir # CHANGE ACCORDINGLY IF ALREADY HAVE THIS DATA
    olr,time,lat,lon,dtime = mync.openolr(dsrc+'olr.day.mean.nc','olr',subs=sub)
    olr=olr[1675:,:,:];time=time[1675:];dtime=dtime[1675:]     
    if testyear: olr, dtime, time = olr[:365,:,:],dtime[:365],time[:365]
    if getdistr: showme = blb.gethists(olr,time,lat,lon,v,sub=sub,figd=figdir)
    # CALLS TO METBLOBS
    plt.ion()
    if getmbs:
        mbs, mbt, chull = blb.MetBlobs(olr,dtime,time,lat,lon,v,\
                                       sub=sub,showblobs=showblb,interact=False)
        blb.mbsave(picdir+v+".mbs",mbs,mbt,chull)
        del mbs,mbt,chull
        if noaaolrall:
            refmbsstr="noaa-olr-0-0"
            refmbs,refmbt,refch = blb.mbopen(picdir+refmbsstr+".mbs")
            reftime=refmbs[:,0]
            v="noaa-olr-0-all"
            dset,varstr, lev, drv = v.split('-')
            exec("ixt,[time,%s,dtime]=\
                  my.ixtwindow(reftime,time,hrwindow,time,%s,dtime)"\
                   %(varstr,varstr) )
            mbs, mbt, chull = blb.MetBlobs(olr,dtime,time,lat,lon,v,\
                                      sub=sub,showblobs=showblb,interact=False)
            blb.mbsave(picdir+v+".mbs",mbs,mbt,chull)
            del mbs,mbt,chull
        if noaaolrfull:
            v="noaa-olr-0-full"
            dset,varstr, lev, drv = v.split('-')
            mbs, mbt, chull = blb.MetBlobs(olr,dtime,time,lat,lon,v,\
                                      sub=sub,showblobs=showblb,interact=False)
            blb.mbsave(picdir+v+".mbs",mbs,mbt,chull)
            del mbs,mbt,chull
        
if synoptics:
    refmbstr="noaa-olr-0-0"
    refall="noaa-olr-0-all"
    reffull="noaa-olr-0-full"
    metblobslist=[refmbstr,refall,reffull]
    mfilelist=[picdir+j+'.mbs' for j in metblobslist]
    s = sy.SynopticEvents(metblobslist,mfilelist,hrwindow=hrwindow)
    s.buildtracks()
    s.buildevents()
    u = s.uniqueevents()
    s.save(picdir+'NOAA-OLR.synop')
    del s

print 'TOTAL TIME TAKEN FOR test.py is:',(tmr.time()-tstart)/60,'mins'
