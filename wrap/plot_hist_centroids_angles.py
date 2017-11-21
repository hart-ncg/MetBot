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


### Running options
cenlonplot=False
cenlatplot=False
angleplot=False
monplot=True

scatter_lon_angle=False
scatter_lat_angle=False
scatter_lon_lat=False

scatter_lon_mon=True
scatter_lat_mon=True
scatter_ang_mon=True

title=True

testyear=True  # plot based on 1 year of test data
testfile=True
threshtest=False # Option to run on thresholds + and - 5Wm2 as a test
cmip5_spec=False

refdset="noaa"
refmod="cdr"

### Directory
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"cen_ang_figs/"

### Loop threshs
if threshtest:
    thnames=['lower','actual','upper']
else:
    thnames=['actual']

nthresh=len(thnames)
for t in range(nthresh):

    ### Multi dset?
    dsets='spec'     # "all" or "spec" to choose specific dset(s)
    #  all doesn't really work because of TRMM
    if dsets=='all':
        ndset=len(dsetdict.dset_deets)
        dsetnames=list(dsetdict.dset_deets)
        dsetstr='all_dset'+'_'+str(ndset)
    elif dsets=='spec': # edit for the dset you want
        #ndset=5
	    #dsetnames=['noaa','ncep','era','20cr','um']
        ndset=1
        dsetnames=['noaa']
        #ndset=6
	    #dsetnames=['noaa','ncep','era','20cr','um','cmip5']
        if cmip5_spec:
            dsetstr = ('_'.join(dsetnames)) + '_spec'
        else:
            dsetstr=('_'.join(dsetnames))+'_'+str(ndset)
    print 'Running on datasets:'
    print dsetnames

    ### Count total number of models
    nm_dset=np.zeros(ndset)
    for d in range(ndset):
        dset = dsetnames[d]
        if cmip5_spec:
            if dset=='cmip5':
                nmod=4
            else:
                nmod = len(dsetdict.dset_deets[dset])
        else:
            nmod = len(dsetdict.dset_deets[dset])
        nm_dset[d]=nmod
    nallmod=np.sum(nm_dset)
    nallmod=int(nallmod)
    print 'Total number of models = '+str(nallmod)

    ### Open array for modnames for cbar
    modnm=["" for x in range(nallmod)] # creates a list of strings for modnames

    ### Display options for plot
    styls = ['solid', 'dashed', 'dotted', 'dashed', 'solid', 'dotted']
    lws = [3, 2, 2, 2, 1, 1]
    zorders = [3, 2, 2, 2, 1, 1]
    mkrs = ["o","^","<",">","x"]
    msiz = [5,5,5,5,5,5]

    ### Open figures
    if cenlonplot: plt.figure(num='cenlon',figsize=[12,10])
    if cenlatplot: plt.figure(num='cenlat',figsize=[12,10])
    if angleplot: plt.figure(num='angle',figsize=[12,10])
    if monplot: plt.figure(num='month',figsize=[12,10])

    if scatter_lon_angle: plt.figure(num='lon_ang',figsize=[12,10])
    if scatter_lat_angle: plt.figure(num='lat_ang', figsize=[12, 10])
    if scatter_lon_lat: plt.figure(num='lon_lat', figsize=[12, 10])

    if scatter_lon_mon: plt.figure(num='lon_mon', figsize=[12, 10])
    if scatter_lat_mon: plt.figure(num='lat_mon', figsize=[12, 10])
    if scatter_ang_mon: plt.figure(num='ang_mon', figsize=[12, 10])

    ### Loop dsets and models
    z=0
    cols = ['k','r', 'b', 'c', 'm', 'g', 'y',\
            'k','r', 'b', 'c', 'm', 'g', 'y',\
            'k','r', 'b', 'c', 'm', 'g', 'y',\
            'k','r', 'b', 'c', 'm', 'g', 'y',\
            'k','r', 'b', 'c', 'm', 'g', 'y']
    for d in range(ndset):
        dset=dsetnames[d]
        if cmip5_spec:
            if dset=='cmip5':
                nmod=4
                mnames=['ACCESS1-0','CanESM2','GFDL-CM3','MIROC-ESM']
            else:
                nmod = len(dsetdict.dset_deets[dset])
                mnames = list(dsetdict.dset_deets[dset])
        else:
            nmod=len(dsetdict.dset_deets[dset])
            mnames=list(dsetdict.dset_deets[dset])
        print 'Looping through models'
        print mnames
        nmstr=str(nmod)

        #	for m in range(1):
        for m in range(nmod):
            name=mnames[m]
            mcnt = str(m + 1)

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

            ### Put name into string list
            modnm[z]=dset+"_"+name

            if cenlonplot:
                plt.figure(num='cenlon')
                y, binEdges = np.histogram(cXs, bins=25, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

            if cenlatplot:
                plt.figure(num='cenlat')
                y, binEdges = np.histogram(cYs, bins=25, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

            if angleplot:
                plt.figure(num='angle')
                y, binEdges = np.histogram(degs, bins=25, density=True)
                bincentres = 0.5 * (binEdges[1:] + binEdges[:-1])
                plt.plot(bincentres, y, linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

            if monplot:
                plt.figure(num='mon')
                season = [8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7]
                mon_count=np.zeros(len(season))
                for imn in xrange(len(mon_count)):
                    mn=season[imn]
                    ix = np.where(mons==mn)[0]
                    mon_count[imn]=len(ix)
                seaspos=np.arange(1,13,1)
                plt.plot(seaspos,mon_count,linestyle=styls[d], linewidth=lws[d], zorder=zorders[d])

            if scatter_lon_angle:
                plt.figure(num='lon_ang')
                plt.scatter(cXs,degs,c=cols[z],marker=mkrs[d],s=msiz[d],edgecolors='face',zorder=zorders[d])

            if scatter_lat_angle:
                plt.figure(num='lat_ang')
                plt.scatter(degs,cYs,c=cols[z],marker=mkrs[d],s=msiz[d],edgecolors='face',zorder=zorders[d])

            if scatter_lon_lat:
                plt.figure(num='lon_lat')
                plt.scatter(cXs,cYs,c=cols[z],marker=mkrs[d],s=msiz[d],edgecolors='face',zorder=zorders[d])

            mon4scatter=np.zeros(len(mons))
            for mn in range(len(mon4scatter)):
                if mons[mn] <=7:
                    mon4scatter[mn]=mons[mn]+5
                elif mons[mn] >=8:
                    mon4scatter[mn]=mons[mn]-7

            if scatter_lon_mon:
                plt.figure(num='lon_mon')
                plt.scatter(cXs,mon4scatter,c=cols[z],marker=mkrs[d],s=msiz[d],edgecolors='face',zorder=zorders[d])

            if scatter_lat_mon:
                plt.figure(num='lat_mon')
                plt.scatter(mon4scatter,cYs,c=cols[z],marker=mkrs[d],s=msiz[d],edgecolors='face',zorder=zorders[d])

            if scatter_ang_mon:
                plt.figure(num='ang_mon')
                plt.scatter(degs,mon4scatter,c=cols[z],marker=mkrs[d],s=msiz[d],edgecolors='face',zorder=zorders[d])

            z += 1


            print 'Finished running on ' + name
            print 'This is model '+mcnt+' of '+nmstr+' in list'

    ### Edits to figures
    if cenlonplot:
        plt.figure(num='cenlon')
        plt.xlim(7.5,100.0)
        plt.xticks([10,20,30,40,50,60,70,80,90,100])
        plt.legend(modnm,loc='upper left',fontsize='xx-small')
        plt.xlabel('Longitude', fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('frequency density', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Distribution of centroid longitudes: '+dsetstr,\
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        cenlonfig=figdir+'/hist_cen_lon.'+dsetstr+'.'+thnames[t]+'.png'
        plt.savefig(cenlonfig)

    if cenlatplot:
        plt.figure(num='cenlat')
        plt.xlim(-15.0,-40.0)
        plt.legend(modnm, loc='upper left', fontsize='xx-small')
        plt.xlabel('Latitude', fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('frequency density', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Distribution of centroid latitudes: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        cenlatfig = figdir + '/hist_cen_lat.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(cenlatfig)

    if angleplot:
        plt.figure(num='angle')
        plt.xlim(-90.0,-5.0)
        plt.legend(modnm, loc='upper left', fontsize='xx-small')
        plt.xlabel('Cloudband Orientation', fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('frequency density', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Distribution of cloudband angles: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        anglefig = figdir + '/hist_angle.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(anglefig)

    monthstr = ['Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']

    if monplot:
        plt.figure(num='mon')
        plt.xticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
        plt.xlim(1,12)
        plt.legend(modnm, loc='upper left', fontsize='xx-small')
        plt.xlabel('Months', fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('frequency', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Total count of CBs by month: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        monfig = figdir + '/hist_mon.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(monfig)


    if scatter_lon_angle:
        plt.figure(num='lon_ang')
        plt.xlim(7.5,100.0)
        plt.ylim(-90.0,-5.0)
        plt.legend(modnm, loc='lower right', fontsize='xx-small')
        plt.xlabel('Longitude of Centroid',fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('Cloudband Orientation', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Relationship between CB longitude and orientation: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        lonangfig = figdir + '/scatter_lon_angle.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(lonangfig)

    if scatter_lat_angle:
        plt.figure(num='lat_ang')
        plt.xlim(-90.0,-5.0)
        plt.ylim(-40.0,-15.0)
        plt.legend(modnm, loc='lower right', fontsize='xx-small')
        plt.xlabel('Cloudband Orientation',fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('Latitude of Centroid', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Relationship between CB orientation and latitude: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        latangfig = figdir + '/scatter_lat_angle.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(latangfig)

    if scatter_lon_lat:
        plt.figure(num='lon_lat')
        plt.xlim(7.5,100.0)
        plt.ylim(-40.0,-15.0)
        plt.legend(modnm, loc='lower right', fontsize='xx-small')
        plt.xlabel('Longitude of Centroid',fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('Latitude of Centroid', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Relationship between CB longitude and latitude: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        lonlatfig = figdir + '/scatter_lon_lat.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(lonlatfig)

    if scatter_lon_mon:
        plt.figure(num='lon_mon')
        plt.xlim(7.5,100.0)
        plt.yticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
        #plt.ylim()
        plt.legend(modnm, loc='lower right', fontsize='xx-small')
        plt.xlabel('Longitude of Centroid',fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('Months', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Relationship between CB longitude and month: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        lonmonfig = figdir + '/scatter_lon_mon.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(lonmonfig)

    if scatter_lat_mon:
        plt.figure(num='lat_mon')
        #plt.xlim()
        plt.xticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
        plt.ylim(-40.0,-15.0)
        plt.legend(modnm, loc='lower right', fontsize='xx-small')
        plt.xlabel('Months',fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('Latitude of Centroid', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Relationship between CB latitude and month: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        latmonfig = figdir + '/scatter_lat_mon.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(latmonfig)

    if scatter_ang_mon:
        plt.figure(num='ang_mon')
        plt.xlim(-90.0,-5.0)
        plt.yticks(np.arange(1, 13), monthstr, fontsize=13.0)  # month labels
        #plt.ylim()
        plt.legend(modnm, loc='lower right', fontsize='xx-small')
        plt.xlabel('Cloudband Orientation',fontsize=14.0, weight='demibold', color='k')
        plt.ylabel('Months', fontsize=14.0, weight='demibold', color='k')
        if title: plt.title('Relationship between CB orientation and month: ' + dsetstr, \
                            fontsize=13.0, weight='demibold', color='k')

        ### Save figure
        angmonfig = figdir + '/scatter_ang_mon.' + dsetstr + '.'+thnames[t]+'.png'
        plt.savefig(angmonfig)


    plt.close('all')

    print 'Finished running on ' + thre_str
    print 'This is thresh ' + str(t) + ' of ' + str(nthresh) + ' in list'
