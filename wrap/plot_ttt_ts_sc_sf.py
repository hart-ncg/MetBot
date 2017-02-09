import iris
iris.FUTURE.netcdf_promote=True
iris.FUTURE.cell_datetime_objects=True
import numpy as np
import matplotlib.pyplot as plt
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
dset="noaa" # options: noaa, um, cmip5
name="noaa" # options: noaa, mo runid, cmip5 model name
ys="1979_1979" # these are the years in the file name
beginatyr="1975" # choose first year for analysis (should fit with metbot run)
vname="olr" # will be olr for most dsets but rlut for cmip5


### Location for olr input & outputs
indir=cwd+"/../../../CTdata/metbot_multi_dset/"+dset+"/"
infile=indir+name+".olr.day.mean."+ys+".nc"
print infile
outdir=indir+name+"/"
outsuf=outdir+name+'_'
mbsfile=outsuf+dset+"-olr-0-0.mbs"
syfile=outsuf+dset+'-OLR.synop'


### Running options
sub="SA"
seasopt="coreseason"    # for spatiofreq plots
                        # options: coreseason, dryseason, fullseason
tsplot=True             # to get timeseries plot
scplot=True             # to get seasonal cycle plots
sfplot=True             # to get spatiofrequency plot
testfile=True           # Uses a test file with short period


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
plt.figure(figsize=[11,5])
plt.plot(yrs,nNF,'k',lw=2.)
plt.savefig(outsuf+dset+'timeseries.png',dpi=150)

### PLOT SEASONAL CYCLE WITH BOX AND WHISKERS
stats.plotseasonbox_rj(scycle,'All__'+count_all,outsuf+dset+'_All',savefig=True)
stats.plotseasonbox_rj(scyclew,'Continental__'+count_cont,outsuf+dset+'_Continental',savefig=True)
stats.plotseasonbox_rj(scyclee,'Oceanic__'+count_mada,outsuf+dset+'_Oceanic',savefig=True)

### PLOT MONTHLY GRIDPOINT COUNT
mapnm=outsuf+dset+'_'+seasopt
msklist=ap.spatiofreq2_season(s,lat,lon,yrs,ks,\
    figno=1,season=seasopt,key=dset+'-olr-0-0',flagonly=True,file_suffix=mapnm,savefig=True)