# MetBlobs wrapper
#
# Designed to be flexible to dataset
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....years of file: $firstyear_lastyear
# .....start year for calculations: currently using 1979 for OLR, 1970 for cmip5
# .....   [adjust depending on years available]
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# Format to follow for all variables
# DSET-VAR-LEV-DERIVE
# DSET-VAR-LEV-DERIVE-{EXP}{ENS} (for flavours experiments [maybe ccam ouput too])
# When time subsets of dsets are used, this should be denoted
# DSET-VAR-LEV-DERIVE-{subsetdescription}
import numpy as np
import datetime
from datetime import date
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
import MetBot.dset_dict as dsetdict
tstart=tmr.time()

### Choose dset and years
dset="cmip5" # options: noaa, um, cmip5
name="ACCESS1-0" # options: noaa, mo runid, cmip5 model name
ys="1975_1975" # these are the years in the file name
beginatyr="1975" # choose first year for analysis
vname="rlut" # will be olr for most dsets but rlut for cmip5

### Location for olr input & outputs
indir=cwd+"/../../../CTdata/metbot_multi_dset/"+dset+"/"
infile=indir+name+".olr.day.mean."+ys+".nc"
print infile
outdir=indir+name+"/"
my.mkdir_p(outdir)
outsuf=outdir+name+'_'

### Running options
olr=True         # Get mbs for $dset-olr-0-0
olrall=True      # Get mbs for $dset-olr-0-all
olrfull=True     # Get mbs for $dset-olr-0-full
testfile=True    # Uses a test file with short period
testyear=True    # Only uses first 365 days of olr data
getdistr=True    # Save a figure showing histogram of OLR values
getmbs=True      # Actually run the MetBot algorithm
showblb=False    # Show the blobs while running
intract=False    # Interactive running of showblobs
refsubset=True   # This is used if noaaolr=True to only look in time window
hrwindow=49      # ... close (49 hours/ 2days) to flagged cloud band days
synoptics=True   # Build tracks of cloud blobs that become TTT cloud bands
                 # ... which are then used to build TTT events.

### Ensure only look at Southern Africa
sub="SA"

### Open OLR data
if olr:
    v=dset+"-olr-0-0"
    daset, varstr, lev, drv = v.split('-')
    ncout = mync.openolr_multi(infile,vname,name,\
                                                dataset=dset,subs=sub)
    ndim = len(ncout)
    if ndim==5:
        olr,time,lat,lon,dtime = ncout
    elif ndim==6:
        olr, time, lat, lon, lev, dtime = ncout
        olr=np.squeeze(olr)
    else:
        print 'Check number of levels in ncfile'

    ### Select data to run
    ### Get time information
    moddct = dsetdict.dset_deets[dset][name]
    units = moddct['timeunit']
    cal = moddct['calendar']
    ### If testfile run on all days available
    if testfile:
        olr = olr[:, :, :];time = time[:];dtime = dtime[:]
    else:
        ### Find starting timestep
        start = moddct['startdate']
        ystart=int(start[0:4]);mstart=int(start[5:7]);dstart=int(start[8:10])
        if cal=="360_day":
            startday=(ystart*360)+((mstart-1)*30)+dstart
            beginday=((int(beginatyr))*360)+1
            daysgap=beginday-startday+1
        else:
            startd=date(ystart,mstart,dstart)
            begind=date(int(beginatyr),01,01)
            daysgap=(begind-startd).days
        olr=olr[daysgap:,:,:];time=time[daysgap:];dtime=dtime[daysgap:]
    if testyear:
        if cal=="360_day":
            olr, dtime, time = olr[:360, :, :], dtime[:360], time[:360]
        else:
            olr, dtime, time = olr[:365,:,:],dtime[:365],time[:365]

    ### Plot olr dist to check threshold
    if getdistr: showme = blb.gethists(olr,time,lat,lon,v,sub=sub,figd=outsuf)
    plt.ion()

    ### Get mbs 0-0
    if getmbs:
        mbs, mbt, chull = blb.MetBlobs(olr,dtime,time,lat,lon,v,\
                                       sub=sub,showblobs=showblb,interact=intract)
        blb.mbsave(outsuf+v+".mbs",mbs,mbt,chull)
        del mbs,mbt,chull

        ### Get mbs 0-all
        if olrall:
            refmbsstr=dset+"-olr-0-0"
            refmbs,refmbt,refch = blb.mbopen(outsuf+refmbsstr+".mbs")
            reftime=refmbs[:,0]
            v=dset+"-olr-0-all"
            daset,varstr, lev, drv = v.split('-')
            exec("ixt,[time,%s,dtime]=\
                  my.ixtwindow(reftime,time,hrwindow,time,%s,dtime)"\
                   %(varstr,varstr) )
            mbs, mbt, chull = blb.MetBlobs(olr,dtime,time,lat,lon,v,\
                                      sub=sub,showblobs=showblb,interact=False)
            blb.mbsave(outsuf+v+".mbs",mbs,mbt,chull)
            del mbs,mbt,chull

        ### Get mbs 0-full
        if olrfull:
            v=dset+"-olr-0-full"
            daset,varstr, lev, drv = v.split('-')
            mbs, mbt, chull = blb.MetBlobs(olr,dtime,time,lat,lon,v,\
                                      sub=sub,showblobs=showblb,interact=False)
            blb.mbsave(outsuf+v+".mbs",mbs,mbt,chull)
            del mbs,mbt,chull

### Get synop file
if synoptics:
    refmbstr=dset+"-olr-0-0"
    refall=dset+"-olr-0-all"
    reffull=dset+"-olr-0-full"
    metblobslist=[refmbstr,refall,reffull]
    mfilelist=[outsuf+j+'.mbs' for j in metblobslist]
    print mfilelist
    s = sy.SynopticEvents(metblobslist,mfilelist,hrwindow=hrwindow)
    s.buildtracks()
    s.buildevents(basetrkkey=refmbsstr)
    u = s.uniqueevents()
    s.save(outsuf+dset+'-OLR.synop')
    del s

print 'TOTAL TIME TAKEN FOR test.py is:',(tmr.time()-tstart)/60,'mins'
