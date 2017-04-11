# Plotting wrapper for OLR comparisons
# options for output
# .... olr histograms (as lines) for multimodels
# .... olr thresholds (from automated fmin detection)
#           as plot
#           as text file
# .... shifted olr distributions by threshold (for illustration)
#
# Option for testfile but will only work if all the files are available
# (use spec to run on certain models only)
#
#  Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5, ncep, era, 20cr
# .....name: noaa or cdr, $mo_runid (e.g. anqjn), $cmip5_model_name, $reanal_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs

### Running options
sub="SA"
seasopt="coreseason"    # for spatiofreq plots
                        # options: coreseason, dryseason, fullseason
histplot=True           # to get olr histograms
threshplot=True         # to get olr threshold plot
threshtext=True         # to put olr thresholds in text file
shiftdist=True          # to plot shifted distributions
testyear=True           # To use output from a test
testfile=True           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
title=True      # plot title
refdset="noaa"
refmod="cdr"
globv='olr'
bkdir=cwd+"/../../../CTdata/metbot_multi_dset"

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['um']
    dsetstr = '_'.join(dsetnames)
ndstr=str(ndset)
print 'Running on datasets:'
print dsetnames

### Get saddle of distribution for ref dset (specified above)
moddct = dsetdict.dset_deets[refdset][refmod]
vname = moddct['olrname']
if testfile:
    ys = moddct['testfileyr']
else:
    ys = moddct['yrfname']
if testyear:
    beginatyr = moddct['startyr']
else:
    beginatyr = moddct['testyr']

indir = bkdir + "/"+ refdset + "/"
infile = indir + refmod + ".olr.day.mean." + ys + ".nc"
print infile
ncout = mync.open_multi(infile, globv, refmod,\
                           dataset=refdset, subs=sub)
ndim = len(ncout)
if ndim == 5:
    olr, time, lat, lon, dtime = ncout
elif ndim == 6:
    olr, time, lat, lon, lev, dtime = ncout
    olr = np.squeeze(olr)
else:
    print 'Check number of levels in ncfile'

### Select data to run
### Get time information
units = moddct['olrtimeunit']
cal = moddct['calendar']
### If testfile run on all days available
if testfile:
    olr = olr[:, :, :];
    time = time[:];
    dtime = dtime[:]
else:
    ### Find starting timestep
    start = moddct['startdate']
    ystart = int(start[0:4]);
    mstart = int(start[5:7]);
    dstart = int(start[8:10])
    if cal == "360_day":
        startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
        beginday = ((int(beginatyr)) * 360) + 1
        daysgap = beginday - startday + 1
    else:
        startd = date(ystart, mstart, dstart)
        begind = date(int(beginatyr), 01, 01)
        daysgap = (begind - startd).days
    olr = olr[daysgap:, :, :]
    time = time[daysgap:]
    dtime = dtime[daysgap:]
if testyear:
    if cal == "360_day":
        olr, dtime, time = olr[:360, :, :], dtime[:360], time[:360]
    else:
        olr, dtime, time = olr[:365, :, :], dtime[:365], time[:365]

### Get saddle of ref olr
refolrvals = olr
refolrthresh=fs.find_saddle(refolrvals,method='fmin',showplot=False)

### Count total number of models - (assumes using "all" models)
nm_dset=np.zeros(ndset)
for d in range(ndset):
    dset = dsetnames[d]
    nmod = len(dsetdict.dset_deets[dset])
    nm_dset[d]=nmod
nallmod=np.sum(nm_dset)
nallmod=int(nallmod)
print 'Total number of models = '+str(nallmod)

### Open array for names for cbar
modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

### Display options for plot
styls=['solid','dashed','dotted','dashed','solid']
lws=[3,2,2,2,1]
zorders=[3,2,2,2,1]

### Open figures
if histplot: plt.figure(num='raw')
if threshtext:
    if testyear:
        txtfile = open(bkdir + "/thresholds.fmin." + dsetstr + ".test.txt", "w")
    else:
        txtfile = open(bkdir + "/thresholds.fmin."+dsetstr+".txt", "w")
if threshplot: plt.figure(num='threshs',figsize=[10,3])
if shiftdist: plt.figure(num='shift')

z=0
### Loop datasets
for d in range(ndset):
    dset=dsetnames[d]
    dcnt=str(d+1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    ### Multi model?
    mods='all'  # "all" or "spec" to choose specific model(s)
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
            beginatyr=moddct['testyr']
        else:
            beginatyr = moddct['startyr']

        ### Location for olr input & outputs
        indir = bkdir + "/" + dset + "/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile

        ### Open olr nc file
        ncout = mync.open_multi(infile,globv,name,\
                                dataset=dset,subs=sub)
        ndim = len(ncout)
        if ndim == 5:
            olr, time, lat, lon, dtime = ncout
        elif ndim == 6:
            olr, time, lat, lon, lev, dtime = ncout
            olr = np.squeeze(olr)
        else:
            print 'Check number of dims in ncfile'

        ### Select olr data
        ### Get time information
        moddct = dsetdict.dset_deets[dset][name]
        units = moddct['olrtimeunit']
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

        ### Get thresh
        olrvals = olr
        olrthresh = fs.find_saddle(olrvals, method='fmin', showplot=False)

        ### Plot histogram with the thresh
        if histplot:
            plt.figure(num='raw')
            olr_flat = np.nan_to_num(olrvals.ravel())
            y, binEdges = np.histogram(olr_flat, bins=50, density=True)
            bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
            plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

        ### Thresh text file
        if threshtext:
            txtfile.write(dset+ "\t" +name+ "\t" + str(int(olrthresh)) + "\n")

        ### Plot thresh
        if threshplot:
            plt.figure(num='threshs')
            plt.plot(olrthresh,1,'^',markersize=20)

        if shiftdist:
            ### Get shifted values
            threshdiff = refolrthresh - olrthresh
            shifted_dist = olr_flat + threshdiff

            ### Plot shifted dist
            plt.figure(num='shift')
            y, binEdges = np.histogram(shifted_dist, bins=50, density=True)
            bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
            plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

        ### Put name into string list
        modnm[z] = dset + "_" + name

        z += 1

        print 'Finished running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'


### Edits to raw figure
if histplot:
    plt.figure(num='raw')

    ### Plot legend and axis
    plt.xlim(100, 320)
    plt.yticks(np.arange(0.002, 0.016, 0.004))
    plt.legend(modnm, loc='upper left',fontsize='xx-small')
    plt.xlabel('OLR', fontsize=13.0, weight='demibold', color='k')
    plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')
    if title: plt.title('Histogram of OLR: '+dsetstr,\
                        fontsize=13.0, weight='demibold', color='k')

    ### Save figure
    rawfig=bkdir+'/olr_raw_hist.'+dsetstr+'.png'
    plt.savefig(rawfig)


### Edits to text file
if threshtext:
    txtfile.close()


### Edits to threshs figure
if threshplot:
    plt.figure(num='threshs')

    ### Add refmean
    plt.plot(refolrthresh,1,'o',c='k',markersize=35,zorder=1)
    plt.xlim(220,260)
    plt.ylim(0,2)
    plt.yticks([0,1,2])
    plt.legend(modnm, loc='upper left',fontsize='xx-small',markerscale=0.5)

    ### Save
    threshfig=bkdir+'/olr_threshs.'+dsetstr+'.png'
    plt.savefig(threshfig)

### Edits to figure with shifted data
if shiftdist:
    plt.figure(num='shift')

    ### Add ref threshold
    plt.plot((refolrthresh, refolrthresh), (0, 0.014),'k')

    ### Plot legend and axes
    plt.xlim(100, 320)
    plt.yticks(np.arange(0.002, 0.016, 0.004))
    plt.legend(modnm, loc='upper left',fontsize='xx-small')
    plt.xlabel('OLR - shifted using '+refdset+'_'+refmod, fontsize=13.0, weight='demibold', color='k')
    plt.ylabel('frequency density', fontsize=13.0, weight='demibold', color='k')
    if title: plt.title('Histogram of shifted OLR: '+dsetstr,\
                        fontsize=13.0, weight='demibold', color='k')

    ### Save figure
    shiffig=bkdir+'/olr_shifted.'+dsetstr+'.png'
    plt.savefig(shiffig)

plt.close('all')