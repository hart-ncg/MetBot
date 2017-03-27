# Scatterplot wrapper
#   to plot number of events for each model versus precipitation
#   for all datasets or specific dataset

# .....dset: noaa, um, cmip5
# .....name: model names will be taken from dset dictionary
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/

import numpy as np
from datetime import date
import matplotlib.pyplot as plt
import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.EventStats as stats
import MetBot.SynopticPlot as syp
import MetBot.AdvancedPlots as ap
import MetBot.RainStats as rs
import MetBot.MetBlobs as blb
import MetBot.mytools as my
import MetBot.mynetcdf as mync
import glob, socket, os
import mpl_toolkits.basemap as bm
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs
import scipy


### Running options
testfile=True  # plot based on test file
testyear=True  # plot based on 1 year of test data
                # will only really work on spec dset & model
threshtest=True # Option to run on thresholds + and - 5Wm2 as a test

### Looping options
prdom=['SA_TR','WPR','EPR'] # Domains over which to average precip
cbdom=['All','Continental','Oceanic'] # CBs from which domain? (based on centroids)
#prdom=['SA_TR','WPR','EPR','WCONT','ECONT'] # Domains over which to average precip
#cbdom=['All','Continental','Oceanic','Wcont','Econt'] # CBs from which domain? (based on centroids)
   # domains are handled together - can't mix and match at the moment
seasons=['ann','djf','jja']

dsets=['all'] # spec or all
mods=['all'] # spec or all - can only be spec if 1 spec dset
if dsets==['spec']:
    ndset=1
    dsetnames=['noaa']
    dsetstr = '_'.join(dsetnames)
    if mods == 'spec':  # edit for the models you want
        nmod = 1
        mnames = ['anqjn']
        nallmod = nmod

### Directory
indir=cwd+"/../../../CTdata/metbot_multi_dset/"
txtdir=indir+"play_output/scatter_nTTT_precip/text/"
picdir=indir+"play_output/scatter_nTTT_precip/plot/"

### Some basics before we start looping
my.mkdir_p(txtdir)
my.mkdir_p(picdir)
globp = 'pr'
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']
nthresh=len(thnames)
if dsets==['all']:
    ndset=len(dsetdict.dset_deets)
    dsetnames=list(dsetdict.dset_deets)
    dsetstr = 'all_dset'
ndstr = str(ndset)
print 'Running on datasets:'
print dsetnames
### Count total number of models - if using all
if mods == ['all']:
    nm_dset = np.zeros(ndset)
    for d in range(ndset):
        dset = dsetnames[d]
        nmod = len(dsetdict.dset_deets[dset])
        nm_dset[d] = nmod
    nallmod = np.sum(nm_dset)
    nallmod = int(nallmod)
    print 'Total number of models = ' + str(nallmod)

### Loop plot options and open arrays
# seasons
for ss in range(len(seasons)):

    # domains
    for r in range(len(cbdom)):

        rain_txtfile = open(txtdir + "precip_mean."\
                       +seasons[ss]+\
                       "."+prdom[r]+".txt", "w")

        # thresholds
        for t in range(nthresh):
            ttt_txtfile = open(txtdir + "ttt_mean." \
                                + thnames[t] + "." + seasons[ss] + \
                                "." + cbdom[r] + ".txt", "w")

### Loop dsets and models
z=0
modnm=["" for x in range(nallmod)]
for d in range(ndset):
    dset=dsetnames[d]
    dcnt = str(d + 1)
    print 'Running on '+dset
    print 'This is dset '+dcnt+' of '+ndstr+' in list'

    if mods == ['all']:
        nmod = len(dsetdict.dset_deets[dset])
        mnames = list(dsetdict.dset_deets[dset])
        nmstr = str(nmod)

    for m in range(nmod):
        name=mnames[m]
        mcnt=str(m+1)
        print 'Running on ' + name
        print 'This is model '+mcnt+' of '+nmstr+' in list'

        ### Open precip file
        if dset == 'noaa':
            raindset = 'trmm'
            rainmod = 'trmm_3b42v7'
        else:
            raindset = dset
            rainmod = name

        raindct = dsetdict.dset_deets[raindset][rainmod]
        rainunits = raindct['prtimeunit']
        raincal = raindct['calendar']
        if testfile:
            rainys = raindct['testfileyr']
        else:
            rainys = raindct['yrfname']
        if testyear:
            rainbeginatyr = raindct['testyr']
        else:
            rainbeginatyr = raindct['startyr']
        rainname = raindct['prname']
        rainfile = indir + raindset + "/" + rainmod + ".pr.day.mean." + rainys + ".nc"
        print rainfile

        # Open precip array
        raincollect=np.zeros((len(prdom),len(seasons)),dtype=np.float32)

        # Loop domains
        for r in range(len(prdom)):

            rainout = mync.open_multi(rainfile, globp, rainmod, \
                                      dataset=raindset, subs=prdom[r])

            ndim = len(rainout)
            if ndim == 5:
                rain, time, rainlat, rainlon, raindates = rainout
            elif ndim == 6:
                rain, time, lat, rainlon, rainlev, raindates = rainout
                rain = np.squeeze(rain)
            else:
                print 'Check number of dims in ncfile'

            ### Select data to run
            ### If testfile run on all days available
            if testfile:
                rain = rain[:, :, :];
                time = time[:];
                raindates = raindates[:]
            else:
                ### Find starting timestep
                start = raindct['startdate']
                ystart = int(start[0:4]);
                mstart = int(start[5:7]);
                dstart = int(start[8:10])
                if raincal == "360_day":
                    startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
                    beginday = ((int(rainbeginatyr)) * 360) + 1
                    daysgap = beginday - startday + 1
                else:
                    startd = date(ystart, mstart, dstart)
                    begind = date(int(rainbeginatyr), 01, 01)
                    daysgap = (begind - startd).days
                rain = rain[daysgap:, :, :];
                time = time[daysgap:];
                raindates = raindates[daysgap:]
            if testyear:
                if raincal == "360_day":
                    rain, raindates, time = rain[:360, :, :], raindates[:360], time[:360]
                else:
                    rain, raindates, time = rain[:365, :, :], raindates[:365], time[:365]


            # Loop seasons
            for ss in range(len(seasons)):

                ### Calculate area average
                # N.B. the "rain" data should now have the right years and domain
                # I am just doing a straight gridbox average
                # ...because it's not spanning too many latitudes
                rainmean=stats.seasmean(rain,raindates,seas=seasons[ss])
                rainmean=round(rainmean, 2)
                raincollect[d,m]=rainmean
                rain_txtfile.write(dset + "\t" + name + "\t" + str(rainmean) + "\n")

        ### Find location of synop file
        outdir = indir + dset + "/" + name + "/"
        if testyear:
            outdir = outdir + 'test/'
        else:
            outdir = outdir
        outsuf = outdir + name + '_'

        ### Get thresh
        threshtxt = indir + 'thresholds.fmin.all_dset.txt'
        #threshtxt = indir + 'thresholds.fmin.noaa.txt'
        print threshtxt
        with open(threshtxt) as f:
            for line in f:
                if dset + '\t' + name in line:
                    thresh = line.split()[2]
                    print 'thresh=' + str(thresh)

        thresh = int(thresh)

        for t in range(nthresh):

            if thnames[t]=='actual':
                thisthresh=thresh
            if thnames[t]=='lower':
                thisthresh=thresh - 5
            if thnames[t]=='upper':
                thisthresh=thresh + 5

            thre_str = str(int(thisthresh))

            ###  Open synop file
            syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'
            s = sy.SynopticEvents((),[syfile],COL=False)

            ### Count number of events
            ks = s.events.keys() # all events
            kw, ke = stats.spatialsubset(s,False,cutlon=40.) # splitting tracks west and east of 40E
            #kww, kwe = stats.spatialsubset(kw,False,cutlon=25.) # splitting tracks west and east of 40E

            # Loop domains
            for r in range(len(cbdom)):
                if cbdom[r]=='All': keys=ks
                if cbdom[r]=='Continental': keys=kw
                if cbdom[r]=='Oceanic': keys = ke
                if cbdom[r]=='WCont': keys = kww
                if cbdom[r]=='ECont': keys = kwe

                scycle, cyclestats, yrs = stats.seasonalcycle(s, keys)

                # Loop seasons
                for ss in range(len(seasons)):

                    if seasons[ss]=='ann': tttmean=scycle.mean
                    if seasons[ss]=='djf': tttmean=scycle[:,4:6].mean
                    if seasons[ss]=='jja': tttmean=scycle[:,(0,10,11)].mean

                    ttt_txtfile.write(dset + "\t" + name + "\t" + str(rainmean) + "\n")

        ### Put name into string list
        modnm[z] = dset + "_" + name

        z += 1

rain_txtfile.close()
ttt_txtfile.close()

# seasons
for ss in range(len(seasons)):

    # domains
    for r in range(len(cbdom)):
        rain_txtfile=txtdir + "precip_mean."\
                    +seasons[ss]+"."+prdom[r]+".txt"

        raindata=np.zeros(nallmod,dtype=float32)

        with open(rain_txtfile) as f:
            for line in f:
                raindata[f] = float(line.split()[2])

        # thresholds
        for t in range(nthresh):
            ttt_txtfile=txtdir + "ttt_mean." \
            + thnames[t] + "." + seasons[ss] + \
            "." + cbdom[r] + ".txt"

            tttdata=np.zeros(nallmod,dtype=float32)

            with open(ttt_txtfile) as f:
                for line in f:
                    tttdata[f] = float(line.split()[2])

            #Plot
            plt.figure()
            plt.scatter(tttdata,raindata)
            grad, int, r_value, p_value, std_err = scipy.stats.linregress(tttdata,raindata)
            rsquared=r_value**2
            plt.plot(tttdata,(grad*tttdata + int),'-')
            plt.xlabel('Number of TTT events')
            plt.ylabel('Mean precip')
            plt.title('Relationship between nTTTs and precip:\n'\
                    +'thresh='+ thnames[t] + "." + seasons[ss] + \
                    "." + cbdom[d] + "."+prdom[d]+"\n"+dsetstr+\
                      ".r2_"+str(round(rsquared,2)))
            scatterfig=picdir+'/scatter_nTTT_precip.thresh='\
                       + thnames[t] + "." + seasons[ss] +\
                       "." + cbdom[r] + "."+prdom[r]+".png"
            plt.savefig(scatterfig,dpi=150)