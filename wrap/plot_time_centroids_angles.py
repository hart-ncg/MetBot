# Wrapper to plot time against centroids and angles

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
from dateutil import rrule
from datetime import datetime


### Which datasets and models?
dsets='spec' # spec or all
if dsets=='spec':
    dsetnames = ['noaa']
    # dsetnames=['noaa','ncep','era','20cr','um']
mods='spec' # spec or all
if mods=='spec':
    mnames = ['cdr']

### Running options
timelonplot=True
minlon=7
maxlon=100

timelatplot=True
minlat=-40
maxlat=-15

timeangplot=True
minang=-90
maxang=-5

group_event=True # If want to plot events as lines True, if individual CBs, false
title=True

testyear=False  # plot based on 1 year of test data
testfile=False
threshtest=False # Option to run on thresholds + and - 5Wm2 as a test
allmodplots=False # Default is to make separate plots for each model,
                    # this option allows one with accumulations from all models to be included

### Directory
bkdir=cwd+"/../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"time_cen_ang_figs/"
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

    # Get time info
    #season = [10, 11, 12, 1, 2, 3, 4]
    #monthstr = ['Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar']
    #nmon = len(season)

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

            ### Count number of flagged CBs
            count_all = len(refmbt)
            print "Total CBs flagged =" + str(count_all)

            ### Open synop file
            if group_event:
                syfile = outsuf + thre_str + '_' + dset + '-OLR.synop'
                s = sy.SynopticEvents((), [syfile], COL=False)
                refkey=dset+'-olr-0-0'
                ks = s.events.keys();ks.sort() # all
                count_all=str(int(len(ks)))
                print "Total CB events ="+str(count_all)


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

            # Loop years
            years=np.unique(edts[:,0])
            print years
            for y in range(len(years)):

                thisyear=years[y]
                if testyear:
                    ix=np.where(edts[:,0]==thisyear)[0]
                else:
                    part1=np.where((edts[:,0]==thisyear) & (edts[:,1]>=8))[0]
                    part2=np.where((edts[:,0]==thisyear+1) & (edts[:,1]<=7))[0]
                    ix=np.concatenate((part1,part2))
                thesedates=edts[ix]

                thesecXs=cXs[ix]
                thesecYs=cYs[ix]
                thesedegs=degs[ix]

                # Find these dates as indices of the full year
                if testyear:
                    firstday=datetime.strptime("01-01-"+str(thisyear),"%d-%m-%Y")
                else:
                    firstday=datetime.strptime("01-08-"+str(thisyear),"%d-%m-%Y")
                fullyear = list(rrule.rrule(rrule.DAILY, count=365, dtstart=firstday))

                full_y_ar=np.zeros([len(fullyear),4],dtype=np.int64)

                for f in range(len(fullyear)):
                    datestr=datetime.strftime(fullyear[f],"%Y-%m-%d-%H")
                    datebits=datestr.split('-')
                    full_y_ar[f,:]=datebits

                # Get labels
                if testyear:
                    monthstr = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul','Aug', 'Sep', 'Oct', 'Nov', 'Dec']
                else:
                    monthstr= ['Aug','Sep','Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar','Apr','May','Jun','Jul']

                monbegins=np.where(full_y_ar[:,2]==1)[0]

                if not group_event:

                    indices=[]
                    for e in range(len(edts)):
                        ix2=my.ixdtimes(full_y_ar,[edts[e][0]],[edts[e][1]],[edts[e][2]],[edts[e][3]])
                        indices.append(ix2)
                    indices=np.squeeze(np.asarray(indices))

                    if timelonplot:

                        plt.figure(figsize=[5,12])
                        plt.scatter(thesecXs,indices)
                        #,marker="o",edgecolour='face',size=20)
                        plt.xlim(7.5,100)
                        plt.ylim(365,0)
                        plt.yticks(monbegins, monthstr, fontsize=13.0)  # month labels
                        plt.xlabel('Centroid Longitude', fontsize=12.0, color='k')
                        if title: plt.title('CB longitudes against time: ' + dset +'_'+name, \
                                        fontsize=12.0, weight='demibold', color='k')

                        timelonfig = figdir + '/time_lon.'+str(thisyear)+'.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                        plt.savefig(timelonfig)
                        plt.close()

                    if timelatplot:

                        plt.figure(figsize=[12,5])
                        plt.scatter(indices,thesecYs)
                        #,marker="o",edgecolour='face',size=20)
                        plt.xlim(0,365)
                        plt.ylim(-40,-15)
                        plt.xticks(monbegins, monthstr, fontsize=13.0)  # month labels
                        plt.xlabel('Centroid Latitude', fontsize=12.0, color='k')
                        if title: plt.title('CB latitudes against time: ' + dset +'_'+name, \
                                        fontsize=12.0, weight='demibold', color='k')

                        timelatfig = figdir + '/time_lat.'+str(thisyear)+'.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                        plt.savefig(timelatfig)
                        plt.close()

                    if timeangplot:

                        plt.figure(figsize=[5,12])
                        plt.scatter(thesedegs,indices)
                        #,marker="o",edgecolour='face',size=20)
                        plt.ylim(365,0)
                        plt.xlim(-90,-5)
                        plt.yticks(monbegins, monthstr, fontsize=13.0)  # month labels
                        plt.xlabel('Cloudband Orientation', fontsize=12.0, color='k')
                        if title: plt.title('CB angles against time: ' + dset +'_'+name, \
                                        fontsize=12.0, weight='demibold', color='k')

                        timeangfig = figdir + '/time_ang.'+str(thisyear)+'.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                        plt.savefig(timeangfig)
                        plt.close()

                elif group_event:

                    if timelonplot:
                        plt.figure(num='time_lon',figsize=[5, 12])

                    if timelatplot:
                        plt.figure(num='time_lat', figsize=[12, 5])

                    if timeangplot:
                        plt.figure(num='time_ang',figsize=[5,12])

                    for k in ks:
                        e = s.events[k]
                        dts = s.blobs[refkey]['mbt'][e.ixflags]
                        indices_ent = []
                        cXs_ent = []
                        cYs_ent = []
                        degs_ent = []

                        for dt in range(len(dts)):
                            x,y = e.trkcX[dt], e.trkcY[dt]
                            cXs_ent.append(x)
                            cYs_ent.append(y)
                            ix3 = my.ixdtimes(full_y_ar, [dts[dt][0]], [dts[dt][1]], [dts[dt][2]], [0])
                            indices_ent.append(ix3)

                            # Find angle
                            ix4 = my.ixdtimes(edts, [dts[dt][0]], [dts[dt][1]], [dts[dt][2]], [0])
                            if len(ix4)==1:
                                degs_ent.append(degs[ix4])
                            elif len(ix4)==0:
                                print 'Date from event not found in mbt...'
                                degs_ent.append(0)
                            elif len(ix4)>1:
                                print 'Date from event found more than once in mbt'
                                all_dates_event=edts[ix4]
                                print all_dates_event
                                all_cXs_event=cXs[ix4]
                                print all_cXs_event
                                print 'this cX='+str(x)
                                ix5=np.where(all_cXs_event==x)[0]
                                right_ind=ix4[ix5]
                                print 'right_ind='+str(right_ind)
                                degs_ent.append(degs[right_ind])

                        indices_ent=np.squeeze(np.asarray(indices_ent))
                        cXs_ent=np.squeeze(np.asarray(cXs_ent))
                        cYs_ent=np.squeeze(np.asarray(cYs_ent))
                        degs_ent=np.squeeze(np.asarray(degs_ent))

                        if timelonplot:
                            plt.figure(num='time_lon')
                            plt.plot(cXs_ent,indices_ent,marker='x')

                        if timelatplot:
                            plt.figure(num='time_lat')
                            plt.plot(indices_ent,cYs_ent,marker='x')

                        if timeangplot:
                            plt.figure(num='time_ang')
                            plt.plot(degs_ent,indices_ent,marker='x')

                    if timelonplot:
                        plt.figure(num='time_lon')
                        plt.xlim(7.5,100)
                        plt.ylim(365,0)
                        plt.xlabel('Centroid Longitude', fontsize=12.0, color='k')
                        plt.yticks(monbegins, monthstr, fontsize=13.0)  # month labels
                        if title: plt.title('CB longitudes against time: ' + dset +'_'+name, \
                                        fontsize=12.0, weight='demibold', color='k')

                        timelonfig = figdir + '/time_lon_events.'+str(thisyear)+'.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                        plt.savefig(timelonfig)
                        plt.close()

                    if timelatplot:
                        plt.figure(num='time_lat')
                        plt.xlim(0,365)
                        plt.ylim(-40,-15)
                        plt.ylabel('Centroid Latitude', fontsize=12.0, color='k')
                        plt.xticks(monbegins, monthstr, fontsize=13.0)  # month labels
                        if title: plt.title('CB latitudes against time: ' + dset +'_'+name, \
                                        fontsize=12.0, weight='demibold', color='k')

                        timelatfig = figdir + '/time_lat_events.'+str(thisyear)+'.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                        plt.savefig(timelatfig)
                        plt.close()

                    if timeangplot:
                        plt.figure(num='time_ang')
                        plt.ylim(365,0)
                        plt.xlim(-90,-5)
                        plt.xlabel('Cloudband Orientation', fontsize=12.0, color='k')
                        plt.yticks(monbegins, monthstr, fontsize=13.0)  # month labels
                        if title: plt.title('CB angles against time: ' + dset +'_'+name, \
                                        fontsize=12.0, weight='demibold', color='k')

                        timeangfig = figdir + '/time_ang_events.'+str(thisyear)+'.' + dset +'_'+name+'.'+ thnames[t] + '.png'
                        plt.savefig(timeangfig)
                        plt.close()


            z += 1

            print 'Finished running on ' + name
            print 'This is model '+mcnt+' of '+nmstr+' in list'


    print 'Finished running on ' + thre_str
    print 'This is thresh ' + str(t) + ' of ' + str(nthresh) + ' in list'
