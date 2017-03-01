# MetBlobs wrapper
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
#
# Option for testfile but will only work if all the files are available
# (use spec to run on certain models only)
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
import os.path
cwd=os.getcwd()
### Add path to MetBot modules and import
sys.path.append(cwd+'/..')
import MetBot.mynetcdf as mync
import MetBot.mytools as my
import MetBot.MetBlobs as blb
import MetBot.SynopticAnatomy as sy
import MetBot.dset_dict as dsetdict
tstart=tmr.time()

### Running options
olr0=True         # Get mbs for $dset-olr-0-0
olrall=True      # Get mbs for $dset-olr-0-all
olrfull=True     # Get mbs for $dset-olr-0-full
testfile=True    # Uses a test file with short period
testyear=True    # Only uses first 365 days of olr data
                 # (testfile designed to be used together with testyear
                 # ..but testyear can be used on any file)
getdistr=True    # Save a figure showing histogram of OLR values
getmbs=True      # Actually run the MetBot algorithm
showblb=False    # Show the blobs while running
intract=False    # Interactive running of showblobs
refsubset=True   # This is used if noaaolr=True to only look in time window
hrwindow=49      # ... close (49 hours/ 2days) to flagged cloud band days
synoptics=True   # Build tracks of cloud blobs that become TTT cloud bands
                 # ... which are then used to build TTT events.
onlynew=False     # Option to only run if the synop file doesn't exist yet

### Ensure only look at Southern Africa
sub="SA"

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['noaa']
ndstr=str(ndset)

for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='spec'  # "all" or "spec" to choose specific model(s)
    if mods=='all':
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
    if mods=='spec': # edit for the models you want
        nmod=1
        mnames=['noaa']
    nmstr=str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        # Get details
        moddct=dsetdict.dset_deets[dset][name]
        vname=moddct['olrname']
        if testfile:
            ys=moddct['testfileyr']
        else:
            ys=moddct['yrfname']
        if testyear:
            beginatyr=moddct['startyr']
        else:
            beginatyr=moddct['testyr']

        ### Location for olr input & outputs
        indir=cwd+"/../../../CTdata/metbot_multi_dset/"+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile
        outdir=indir+name+"/"
        if testyear: outdir=outdir+'test/'
        else: outdir=outdir
        my.mkdir_p(outdir)
        outsuf=outdir+name+'_'

        # Check if file exists for this model
        if onlynew:
            chfile = outsuf + dset + '-OLR.synop'
            if os.path.isfile(chfile):
                print "MetBot already run on this model: "+name
                continue # goes back to the beginning of the for loop
            else:
                print "Running for the first time on: "+name

        ### Open OLR data
        if olr0:
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

        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

print 'TOTAL TIME TAKEN FOR get_ttcbs_loopmodel.py is:',(tmr.time()-tstart)/60,'mins'
