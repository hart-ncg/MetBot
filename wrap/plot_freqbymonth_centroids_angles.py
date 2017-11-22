# Wrapper to plot distribution of centroids and angles from different models
#   for all datasets or specific dataset
#   now also has some scatter plot options to plot relationship between centroid and angle

# .....dset: noaa, ncep, era, 20cr, um, cmip5
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


### Which datasets and models?
dsets='spec' # spec or all
if dsets=='spec':
    dsetnames = ['noaa']
    # dsetnames=['noaa','ncep','era','20cr','um']
mods='spec' # spec or all
if mods=='spec':
    mnames = ['cdr']

### Running options
freqlonplot=True
minlon=7
maxlon=100

freqlatplot=True
minlat=-40
maxlat=-15

freqangplot=True
minang=-90
maxang=-5

title=True

testyear=False  # plot based on 1 year of test data
testfile=False
threshtest=False # Option to run on thresholds + and - 5Wm2 as a test
allmodplots=False # Default is to make separate plots for each model,
                    # this option allows one with accumulations from all models to be included

### Directory
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"freqbymon_figs/"
my.mkdir_p(figdir)

### Loop threshs
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']

nthresh=len(thnames)
for t in range(nthresh):

    if dsets=='all':
        ndset=len(dsetdict.dset_deets)
        dsetnames=list(dsetdict.dset_deets)
        dsetstr='all_dset'+'_'+str(ndset)
    elif dsets=='spec': # edit for the dset you want
        ndset=len(dsetnames)
        dsetstr=('_'.join(dsetnames))+'_'+str(ndset)
    ndstr = str(ndset)
    print 'Running on datasets:'
    print dsetnames

    ### Count total number of models
    nm_dset=np.zeros(ndset)
    for d in range(ndset):
        dcnt = str(d + 1)
        dset = dsetnames[d]
        nmod = len(dsetdict.dset_deets[dset])
        nm_dset[d]=nmod
    nallmod=np.sum(nm_dset)
    nallmod=int(nallmod)
    print 'Total number of models = '+str(nallmod)

    ### Get info for arrays
    nlon = maxlon - minlon
    lons=np.arange(minlon,maxlon+1,1)

    nlat = maxlat - minlat
    lats=np.arange(minlat,maxlat+1,1)

    nang = maxang - minang
    angs=np.arange(minang,maxang+1,1)

    season = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
    season4mesh = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']
    nmon = len(season)

    ### If allmodplot open array for all models
    if allmodplots:
        if freqlonplot:
            all_lon_array=np.zeros([nlon,nmon,nallmod])
        if freqlatplot:
            all_lat_array=np.zeros([nmon,nlat,nallmod])
        if freqangplot:
            all_ang_array=np.zeros([nang,nmon,nallmod])

    ### Loop dsets and models
    z=0
    for d in range(ndset):
        dset=dsetnames[d]
        dcnt=str(d+1)
        print 'Running on '+dset
        print 'This is dset '+dcnt+' of '+ndstr+' in list'

        ### Loop models
        if mods == 'all':
            nmod = len(dsetdict.dset_deets[dset])
            mnames = list(dsetdict.dset_deets[dset])
        if mods == 'spec':  # edit for the models you want
            nmod = len(mnames)
            nmstr = str(nmod)
        print 'Looping through models'
        print mnames
        nmstr=str(nmod)

        for m in range(nmod):
            name=mnames[m]
            mcnt = str(m + 1)

            print 'Running on ' + name
            print 'This is model ' + mcnt + ' of ' + nmstr + ' in list'

            ### Find location of mbs file
            sydir=botdir+dset+"/"+name+"/"
            if testyear: outdir=sydir+'test/'
            else: outdir=sydir
            outsuf=outdir+name+'_'

            ### Get thresh
            if testyear:
                threshtxt = botdir + 'thresholds.fmin.' + dset + '.test.txt'
            else:
                threshtxt = botdir + 'thresholds.fmin.all_dset.txt'
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

            ###  Open mbs file
            mbsfile = outsuf + thre_str + '_' + dset + "-olr-0-0.mbs"
            refmbs, refmbt, refch = blb.mbopen(mbsfile)
            refmbt[:, 3] = 0

            ### Count number of events
            count_all = len(refmbt)
            print "Total CBs flagged =" + str(count_all)

            ### Get array of centroids and angles
            edts = []
            cXs = []
            cYs = []
            degs = []
            mons = []

            for b in range(len(refmbt)):
                date = refmbt[b]
                mon = int(date[1])
                cX = refmbs[b, 3]
                cY = refmbs[b, 4]
                deg = refmbs[b, 2]

                edts.append(date)
                cXs.append(cX)
                cYs.append(cY)
                degs.append(deg)
                mons.append(mon)

            edts = np.asarray(edts)
            edts[:, 3] = 0
            cXs = np.asarray(cXs)
            cYs = np.asarray(cYs)
            degs = np.asarray(degs)
            mons = np.asarray(mons)

            ### Get month by longitude array
            mon_lon_count=np.zeros([nlon,nmon])
            mon_lat_count=np.zeros([nmon,nlat])
            mon_ang_count=np.zeros([nang,nmon])

            for imn in xrange(nmon):
                mn = season[imn]
                ix = np.where(mons == mn)[0]
                cXs_thismon=cXs[ix]
                cYs_thismon=cYs[ix]
                degs_thismon=degs[ix]

                for lo in xrange(nlon):
                    ln=lons[lo]
                    ix2=np.where((cXs_thismon>=ln) & (cXs_thismon <ln+1))[0]

                    mon_lon_count[lo,imn]=len(ix2)
                    if allmodplots:
                        all_lon_array[lo,imn,z]=len(ix2)

                for la in xrange(nlat):
                    lt=lats[la]
                    ix3=np.where((cYs_thismon>=lt) & (cYs_thismon <lt+1))[0]

                    mon_lat_count[imn,la]=len(ix3)
                    if allmodplots:
                        all_lat_array[imn,la,z]=len(ix3)

                for an in xrange(nang):
                    ag=angs[an]
                    ix4=np.where((degs_thismon>=ag) & (degs_thismon<ag+1))[0]

                    mon_ang_count[an,imn]=len(ix4)
                    if allmodplots:
                        all_ang_array[an,imn,z]=len(ix4)

            if freqlonplot:

                plt.figure(figsize=[8,5])
                pmon,plon=np.meshgrid(season4mesh,lons)
                cm=plt.cm.viridis
                cs=plt.pcolormesh(plon,pmon,mon_lon_count,cmap=cm)
                plt.yticks(np.arange(0.5, 12.5, 1), monthstr, fontsize=13.0)  # month labels
                plt.ylim(0, 12)
                plt.xlim(7,100)
                plt.colorbar()
                plt.xlabel('Centroid Longitude', fontsize=12.0, color='k')
                if title: plt.title('Frequency of CB longitudes by month: ' + dsetstr, \
                                    fontsize=12.0, weight='demibold', color='k')

                freqlonfig = figdir + '/freqbymon_lon.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                plt.savefig(freqlonfig)
                plt.close()

            if freqlatplot:

                plt.figure(figsize=[8,5])
                plat,pmon=np.meshgrid(lats,season4mesh)
                cm=plt.cm.viridis
                cs=plt.pcolormesh(pmon,plat,mon_lat_count,cmap=cm)
                plt.xticks(np.arange(0.5, 12.5, 1), monthstr, fontsize=13.0)  # month labels
                plt.xlim(0, 12)
                plt.ylim(-40,-15)
                plt.colorbar()

                plt.ylabel('Centroid Latitude', fontsize=12.0, color='k')
                if title: plt.title('Frequency of CB latitudes by month: ' + dsetstr, \
                                    fontsize=12.0, weight='demibold', color='k')


                freqlatfig = figdir + '/freqbymon_lat.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                plt.savefig(freqlatfig)
                plt.close()

            if freqangplot:

                plt.figure(figsize=[8,5])
                pmon,pang=np.meshgrid(season4mesh,angs)
                cm=plt.cm.viridis
                cs=plt.pcolormesh(pang,pmon,mon_ang_count,cmap=cm)
                plt.yticks(np.arange(0.5, 12.5, 1), monthstr, fontsize=13.0)  # month labels
                plt.ylim(0, 12)
                plt.xlim(-90,-5)
                plt.colorbar()

                plt.xlabel('Centroid Orientation', fontsize=12.0, color='k')
                if title: plt.title('Frequency of CB angles by month: ' + dsetstr, \
                                    fontsize=12.0, weight='demibold', color='k')


                freqangfig = figdir + '/freqbymon_ang.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                plt.savefig(freqangfig)
                plt.close()


            z += 1

            print 'Finished running on ' + name
            print 'This is model '+mcnt+' of '+nmstr+' in list'


    print 'Finished running on ' + thre_str
    print 'This is thresh ' + str(t) + ' of ' + str(nthresh) + ' in list'
