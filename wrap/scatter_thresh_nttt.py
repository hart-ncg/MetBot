# Scatterplot wrapper
#   to plot number of events for each model versus threshold
#   for all datasets or specific dataset

# .....dset: noaa, um, cmip5
# .....name: model names will be taken from dset dictionary
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/

import numpy as np
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
testyear=False  # plot based on 1 year of test data
threshtest=False # Option to run on thresholds + and - 5Wm2 as a test

### Directory
indir=cwd+"/../../../CTdata/metbot_multi_dset/"

### Loop threshs
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']

nthresh=len(thnames)

for t in range(nthresh):

    ### Multi dset?
    dsets='spec'     # "all" or "spec" to choose specific dset(s)
    if dsets=='all':
        ndset=len(dsetdict.dset_deets)
        dsetnames=list(dsetdict.dset_deets)
        dsetstr='all_dset'+'_'+str(ndset)
    elif dsets=='spec': # edit for the dset you want
	ndset=6
	dsetnames=['noaa','ncep','era','20cr','um','cmip5']
#        ndset=2
#        dsetnames=['noaa','um']
        dsetstr=('_'.join(dsetnames))+'_'+str(ndset)
    print 'Running on datasets:'
    print dsetnames


    ### Count total number of models
    nm_dset=np.zeros(ndset)
    for d in range(ndset):
        dset = dsetnames[d]
        nmod = len(dsetdict.dset_deets[dset])
        nm_dset[d]=nmod
    nallmod=np.sum(nm_dset)
    nallmod=int(nallmod)
    print 'Total number of models = '+str(nallmod)

    ### Open arrays for results
    ttt_count=np.zeros((3,nallmod))
    threshlist=np.zeros(nallmod)
    modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

    ### Loop dsets and models
    z=0
    for d in range(ndset):
        dset=dsetnames[d]
        nmod=len(dsetdict.dset_deets[dset])
        mnames=list(dsetdict.dset_deets[dset])
        print 'Looping through models'
        print mnames

        for m in range(nmod):
            name=mnames[m]

            ### Find location of synop file
            outdir=indir+dset+"/"+name+"/"
            if testyear: outdir=outdir+'test/'
            else: outdir=outdir
            outsuf=outdir+name+'_'

            ### Get thresh
            if testyear:
                threshtxt = indir + 'thresholds.fmin.'+dset+'.test.txt'
            else:
                threshtxt=indir+'thresholds.fmin.all_dset.txt'
            print threshtxt
            with open(threshtxt) as f:
                for line in f:
                    if dset+'\t'+name in line:
                        thresh = line.split()[2]
                        print 'thresh='+str(thresh)

            thresh = int(thresh)

            if thnames[t]=='actual':
                thisthresh=thresh
            if thnames[t]=='lower':
                thisthresh=thresh - 5
            if thnames[t]=='upper':
                thisthresh=thresh + 5

            thre_str = str(int(thisthresh))
            threshlist[z] = thisthresh

            ###  Open synop file
            syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'
            s = sy.SynopticEvents((),[syfile],COL=False)

            ### Count number of events
            ks = s.events.keys() # all events
            kw, ke = stats.spatialsubset(s,False,cutlon=40.) # splitting tracks west and east of 40E

            ### Put n events into array
            ttt_count[0,z]=len(ks)
            ttt_count[1,z]=len(kw)
            ttt_count[2,z]=len(ke)

            ### Put name into string list
            modnm[z]=dset+"_"+name
            z+=1

    ### Loop domains
    doms=['All','Continental','Madagascar']
    ndoms=len(doms)

    for r in range(ndoms):
        val=ttt_count[r,:]
        val4plot=val
        mod4plot=modnm

        #for m in range(nallmod):

        plt.figure()
        plt.scatter(threshlist,val4plot)
        m, c, r_value, p_value, std_err = scipy.stats.linregress(threshlist,val4plot)
        rsquared=r_value**2
        plt.plot(threshlist,(m*threshlist + c),'-')
        plt.xlim(230,265)
        plt.xlabel('OLR threshold')
        plt.ylabel('Number of TTT events')
        plt.title(doms[r]+'.thresh_'+thnames[t]+'_r2_'+str(round(rsquared,2)))
        scatterfig=indir+'/scatter_threshold_nTTT.thresh_'+thnames[t]+'.'+dsetstr+'.'+doms[r]+'.png'
        plt.savefig(scatterfig,dpi=150)
