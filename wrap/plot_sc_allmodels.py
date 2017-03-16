# Plotting wrapper
# to plot
# .... box plots for seasonal cycle of TTTs
# ..... to show spread between models
# ....   for whole domain, and west & east of 40E
#
# Comparing:
#       original OLR thresholds
#       OLR thresholds detected automatically using "find_saddle"
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
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
import MetBot.find_saddle as fs
#from MetBot.mytools import savef
import mpl_toolkits.basemap as bm


### Directory
indir=cwd+"/../../../CTdata/metbot_multi_dset/"
picdir=indir+'season_thresh_figs/'


### Running options
sub="SA"
monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']
meanplot=True             # to get boxplot based on mean of years
medplot=True              # to get boxplot based on median of years
allyearplot=True          # to get boxplot based on all years

### Multi dset?
dsets='spec'     # "all" or "spec" to choose specific dset(s)
if dsets=='all':
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset' + '_' + str(ndset)
elif dsets=='spec': # edit for the dset you want
    ndset=1
    dsetnames=['cmip5']
    dsetstr = ('_'.join(dsetnames)) + '_' + str(ndset)
ndstr=str(ndset)

### Count total number of models
nm_dset = np.zeros(ndset)
for d in range(ndset):
    dset = dsetnames[d]
    nmod = len(dsetdict.dset_deets[dset])
    nm_dset[d] = nmod
nallmod = np.sum(nm_dset)
nallmod = int(nallmod)
print 'Total number of models = ' + str(nallmod)

### Loop auto or man threshold
th_type=['man','auto']
for q in range(len(th_type)):

    ### Open arrays for results
    if meanplot:
        mean_data=np.zeros((3,nallmod,12))
    if medplot:
        med_data=np.zeros((3,nallmod,12))
    if allyearplot:
        all_data=np.zeros((3,nallmod*35,12))

    modnm = ["" for x in range(nallmod)]  # creates a list of strings for modnames

    z=0
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

            ### Find location of synop file
            outdir = indir + dset + "/" + name + "/"
            outsuf = outdir + name + '_'

            if th_type[q]=='man':
                ###  Open synop file
                syfile = outsuf + '_' + dset + '-OLR.synop'
                s = sy.SynopticEvents((), [syfile], COL=False)

            elif th_type[q]=='auto':

                ### Get thresh
                with open(indir + 'thresholds.fmin.all_dset.txt') as f:
                    for line in f:
                        if dset + '\t' + name in line:
                            thresh = line.split()[2]
                            print 'thresh=' + str(thresh)

                thisthresh = int(thresh)
                thre_str = str(thisthresh)

                ###  Open synop file
                syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'
                s = sy.SynopticEvents((),[syfile],COL=False)

            ### Count number of events
            ks = s.events.keys();ks.sort()  # all events
            kw, ke = stats.spatialsubset(s, False, cutlon=40.)  # splitting tracks west and east of 40E

            ### Calc seasonal cycle
            scycle, cyclestats, yrs = stats.seasonalcycle(s,False)
            scyclew, cyclestatsw, yrsw = stats.seasonalcycle(s,kw)
            scyclee, cyclestatse, yrse = stats.seasonalcycle(s,ke)

            ### Loop domains
            doms = ['All', 'Continental', 'Madagascar']
            ndoms = len(doms)
            whichc = [sycle,scyclew,scyclee]

            for r in range(ndoms):

                ### Put into data arrays
                if meanplot:
                    monmean=whichc[r].mean()
                    mean_data[r,z,:] = monmean
                if medplot:
                    monmed=whichc[r].median()
                    med_data[r,z,:] = monmed
                if allyearplot:
                    firstdate=z*35
                    lastdate=firstdate+34
                    all_data[r,firstdate:lastdate,:] = whichc[r]

            ### Put name into string list
            modnm[z] = dset + "_" + name
            z += 1

    ### Now we have arrays, time to plot...
    for r in range(ndoms):

        if meanplot:

            plotvals=np.squeeze(mean_data[r,:,:])
            plt.figure()
            plt.boxplot(plotvals, notch=0, sym='+', vert=1, whis=1.5)  # produces boxplot
            plt.plot(np.arange(1, 13), plotvals.mean(0), 'k-', lw=1)  # produces mean of model means
            plt.xticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
            plt.yticks(np.arange(1, 14), fontsize=13.0)
            plt.ylim(0, 8.5)
            plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
            plttitle='Intra-model spread of monthly mean CBs: with '\
                     +th_type[q]+' thresh_'+doms[r]
            plt.title(plttitle, fontweight='demibold')
            # plt.grid()

            fname = picdir+'Intra-model_spread.mean.'\
                    +th_type[q]+'.'+doms[r]+'.'+dsetstr+'.png'
            if savefig: plt.savefig(fname, dpi=150)

        if medplot:

            plotvals=np.squeeze(med_data[r,:,:])
            plt.figure()
            plt.boxplot(plotvals, notch=0, sym='+', vert=1, whis=1.5)  # produces boxplot
            plt.plot(np.arange(1, 13), plotvals.mean(0), 'k-', lw=1)  # produces mean of model medians
            plt.xticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
            plt.yticks(np.arange(1, 14), fontsize=13.0)
            plt.ylim(0, 8.5)
            plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
            plttitle='Intra-model spread of monthly median CBs: with '\
                     +th_type[q]+' thresh_'+doms[r]
            plt.title(plttitle, fontweight='demibold')
            # plt.grid()

            fname = picdir+'Intra-model_spread.median.'\
                    +th_type[q]+'.'+doms[r]+'.'+dsetstr+'.png'
            if savefig: plt.savefig(fname, dpi=150)

        if allyearplot:

            plotvals=np.squeeze(all_data[r,:,:])
            plt.figure()
            plt.boxplot(plotvals, notch=0, sym='+', vert=1, whis=1.5)  # produces boxplot
            plt.plot(np.arange(1, 13), plotvals.mean(0), 'k-', lw=1)  # produces mean of all data
            plt.xticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
            plt.yticks(np.arange(1, 14), fontsize=13.0)
            plt.ylim(0, 8.5)
            plt.ylabel('No. of Cloudbands', fontsize=13.0, weight='demibold')
            plttitle='Intra-model spread of monthly CBs - all years: with '\
                     +th_type[q]+' thresh_'+doms[r]
            plt.title(plttitle, fontweight='demibold')
            # plt.grid()

            fname = picdir+'Intra-model_spread.allyear.'\
                    +th_type[q]+'.'+doms[r]+'.'+dsetstr+'.png'
            if savefig: plt.savefig(fname, dpi=150)