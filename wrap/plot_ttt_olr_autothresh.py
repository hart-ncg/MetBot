# Plotting wrapper
# to plot
# ....olr vals associated with TTTs
# .... aim to explore the relationship with convection in tropics a bit
#
# OLR threshold is detected automatically using "find_saddle"
# Option to run on other OLR thresholds a test - currently + and - 5Wm2
#
#
# Designed to be flexible to dataset
# and run on multiple models in a loop
# input at the top
# .....dset: noaa, um, cmip5
# .....name: noaa, $mo_runid (e.g. anqjn), $cmip5_model_name
# .....directory: here ../../CTdata/metbot_multi_dset/$dset/
# naming of ncfiles used here /$dset/$name.olr.day.mean.$firstyear_$lastyear.nc
# but at present time have focused on trmm and noaa so the model stuff might not work well


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
import mpl_toolkits.basemap as bm

### Running options
sub="SA"
#sub="SA_CONT"
#sub="UM_FOC"
seasopt="fullseason"    # options: coreseason, dryseason, fullseason
testyear=True           # To use output from a test
testfile=True           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
threshtest=False         # Option to run on thresholds + and - 5Wm2 as a test

allplot=False            # plot ave OLR
ave_ttt_plot=False      # plot ave OLR on TTT days - composite
comp_anom_ttt_plot=False  # plot ave OLR on TTT days as anom from long term daily mean for each month
comp_anom_ag_plot=False   # plot comp anom with agtest on composite
comp_anom_cnt_plot=True     # plot count of the number of days above or below average
perc_ag=70              # show if this % or more days agree

under_dayof='dayof'     # if "dayof" plots OLR on TTT days
                        #   if "under" plots olr under TTTs (based on blobs)
monmean='day'           # 'day' is daily mean
                        # 'mon' is monthly mean
nTTTlab=True            # labels each plot with # or % of TTTs

freecol=False            # free colour bar
refkey='0'              # 0 or all
#doms=['All']
doms=['All','Cont','Mada'] # doms for TTT days selected


bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
olrdir=bkdir+"olr_figs/"

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
        mnames=['cdr']
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
        units = moddct['olrtimeunit']
        cal = moddct['calendar']

        ### Location for olr input & outputs
        indir=bkdir+dset+"/"
        infile=indir+name+".olr.day.mean."+ys+".nc"
        print infile
        outdir=indir+name+"/"
        if testyear: outdir=outdir+'test/'
        else: outdir=outdir
        my.mkdir_p(outdir)
        outsuf=outdir+name+'_'

        ### Open olr nc file
        v = dset + "-olr-0-0"
        daset, globv, lev, drv = v.split('-')
        ncout = mync.open_multi(infile,globv,name,\
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

        ### Get thresholds and loop
        if testyear:
            threshtxt = bkdir + 'thresholds.fmin.'+dset+'.test.txt'
        else:
            threshtxt=bkdir+'thresholds.fmin.all_dset.txt'
        print threshtxt
        with open(threshtxt) as f:
            for line in f:
                if dset+'\t'+name in line:
                    thresh = line.split()[2]
                    print 'thresh='+str(thresh)
        thresh = int(thresh)
        if threshtest:
            lowert = thresh - 5
            uppert = thresh + 5
            threshs = [lowert, thresh, uppert]
        else:
            threshs = [thresh]
        ### Loop threshes
        nthresh=len(threshs)
        for t in range(nthresh):
            thisthresh=threshs[t]
            print thisthresh
            thre_str=str(int(thisthresh))
            print thre_str

            mbsfile=outsuf+thre_str+'_'+dset+"-olr-0-0.mbs"
            syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'

            ### Open ttt data
            s = sy.SynopticEvents((),[syfile],COL=False)
            refmbs, refmbt, refch = blb.mbopen(mbsfile)

            ### Select events
            ks = s.events.keys();ks.sort() # all
            count_all=str(int(len(ks)))
            print "Total CB events ="+str(count_all)
            kw, ke = stats.spatialsubset(s,False,cutlon=40.) # events west and east of 40E
            keys=[ks,kw,ke]

            ### Plot olrmaps
            olrbase=olrdir+dset+"/"
            my.mkdir_p(olrbase)
            mapsuf = seasopt+'_'+sub+'_'+dset+'_'+name+'_'+thre_str+'_key'+refkey+'_4'+monmean
            if testfile or testyear:
                testq=True
            else:
                testq=freecol


            if allplot:
                # Only plot this for one threshold and one dom
                if t==0:
                    print 'Plotting ave olr'
                    msklist=ap.gridolrmap_season(s,ks,olr,lat,lon,dtime,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='ave_all',mmean=monmean,\
                                                  under_of=under_dayof,figdir=olrbase,file_suffix=mapsuf,\
                                                  savefig=True,test=testq)

            # Loop domains
            for do in range(len(doms)):
                domname=doms[do]
                print "Running on domain: "+domname
                eventkeys=keys[do]

                newsuf=mapsuf+'_'+domname

                if ave_ttt_plot:
                    print 'Plotting ave olr from TTTs'
                    msklist=ap.gridolrmap_season(s,eventkeys,olr,lat,lon,dtime,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='ave_ttt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=olrbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab)

                if comp_anom_ttt_plot:
                    print 'Plotting ave olr from TTTs'
                    msklist=ap.gridolrmap_season(s,eventkeys,olr,lat,lon,dtime,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ttt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=olrbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab)

                if comp_anom_ag_plot:
                    newnewsuf=newsuf+'_agthr'+str(perc_ag)
                    print 'Plotting composite rainfall anomalies with ag test'
                    msklist=ap.gridolrmap_season(s,eventkeys,olr,lat,lon,dtime,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ag',mmean=monmean,\
                                                  under_of=under_dayof,figdir=olrbase,file_suffix=newnewsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab, agthresh=perc_ag)

                if comp_anom_cnt_plot:
                    msklist=ap.gridolrmap_season(s,eventkeys,olr,lat,lon,dtime,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_cnt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=olrbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq)