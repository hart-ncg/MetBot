# Plotting wrapper
# to plot
# .... time series of no of core season TTTs over time
# .... box plots for seasonal cycle of TTTs
# ....   for whole domain, and west & east of 40E
# .... gridpoint frequency maps for each month
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
import iris
iris.FUTURE.netcdf_promote=True
iris.FUTURE.cell_datetime_objects=True
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.AdvancedPlots as ap
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import MetBot.dset_dict as dsetdict
#from MetBot.mytools import savef
import mpl_toolkits.basemap as bm


### Choose dset and years
dset="cmip5" # options: noaa, um, cmip5
name="ACCESS1-0" # options: noaa, mo runid, cmip5 model name
ys="1975_1975" # these are the years in the file name
beginatyr="1975" # choose first year for analysis (should fit with metbot run)
vname="rlut" # will be olr for most dsets but rlut for cmip5

### Running options
sub="SA"
seasopt="coreseason"    # for spatiofreq plots
                        # options: coreseason, dryseason, fullseason
tsplot=True             # to get timeseries plot
scplot=True             # to get seasonal cycle plots
sfplot=True             # to get spatiofrequency plot
testyear=True           # To use output from a test
testfile=True           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
res='native'            # Option to plot at 'noaa' res or 'native' res

### Location for olr input & outputs
indir=cwd+"/../../../CTdata/metbot_multi_dset/"+dset+"/"
infile=indir+name+".olr.day.mean."+ys+".nc"
print infile
outdir=indir+name+"/"
if testyear:outdir = outdir + 'test/'
else:outdir = outdir
outsuf=outdir+name+'_'
mbsfile=outsuf+dset+"-olr-0-0.mbs"
syfile=outsuf+dset+'-OLR.synop'

### Open olr nc file
ncout = mync.openolr_multi(infile,vname,name,\
                        dataset=dset,subs=sub)
ndim = len(ncout)
if ndim == 5:
    olr, time, lat, lon, dtime = ncout
elif ndim == 6:
    olr, time, lat, lon, lev, dtime = ncout
    olr = np.squeeze(olr)
else:
    print 'Check number of levels in ncfile'


### Select olr data
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

### Option to open noaa file for res
if res=='native':
    lat2=lat
    lon2=lon
elif res=='noaa':
    if testfile:
        yr_noaa="1979_1979"
    else:
        yr_noaa="1974_2013"
    f_noaa=cwd+"/../../../CTdata/metbot_multi_dset/"\
        "noaa/noaa.olr.day.mean."+yr_noaa+".nc"
    olrdump,timedump,lat2,lon2,dtimedump = mync.openolr(f_noaa,'olr',subs=sub)

### Open ttt data
s = sy.SynopticEvents((),[syfile],COL=False)
refmbs, refmbt, refch = blb.mbopen(mbsfile)

### Select events
ks = s.events.keys();ks.sort() # all
kw, ke = stats.spatialsubset(s,False,cutlon=40.) # events west and east of 40E

### Count number of events
count_all=str(int(len(ks)))
count_cont=str(int(len(kw)))
count_mada=str(int(len(ke)))

### Calc seasonal cycle
scycle, cyclestats, yrs = stats.seasonalcycle(s,False)
scyclew, cyclestatsw, yrsw = stats.seasonalcycle(s,kw)
scyclee, cyclestatse, yrse = stats.seasonalcycle(s,ke)
nNF=scycle[:,3:7].sum(1) # summing years for months November to March

### PLOT TIMESERIES OF SEASONAL CYCLE
print 'Plotting timeseries'
plt.figure(figsize=[11,5])
plt.plot(yrs,nNF,'k',lw=2.)
plt.savefig(outsuf+dset+'_timeseries.png',dpi=150)

### PLOT SEASONAL CYCLE WITH BOX AND WHISKERS
print 'Plotting seasonal cycle'
stats.plotseasonbox_rj(scycle,'All__'+count_all,outsuf+dset+'_All',savefig=True)
stats.plotseasonbox_rj(scyclew,'Continental__'+count_cont,outsuf+dset+'_Continental',savefig=True)
stats.plotseasonbox_rj(scyclee,'Oceanic__'+count_mada,outsuf+dset+'_Oceanic',savefig=True)

### PLOT MONTHLY GRIDPOINT COUNT
print 'Plotting spatiofrequency'
mapnm=outsuf+dset+'_'+seasopt
msklist=ap.spatiofreq3_season(s,lat2,lon2,yrs,ks,\
    figno=1,season=seasopt,key=dset+'-olr-0-0',res=res,flagonly=True,file_suffix=mapnm,savefig=True)