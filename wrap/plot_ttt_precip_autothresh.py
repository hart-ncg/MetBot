# Plotting wrapper
# to plot
# ....precip associated with TTTs
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
subrain="SA_TRMM"
#subrain="SA_CONT"
#subrain="UM_FOC"
seasopt="coreseason"    # options: coreseason, dryseason, fullseason
testyear=False           # To use output from a test
testfile=False           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
threshtest=False         # Option to run on thresholds + and - 5Wm2 as a test

# options to plot all days
allplot=False            # plot total rainfall
all_cnt_anom=False       # plot % of days which have +ve anomalies

# options to plot ttt...
tot_ttt_plot=False      # plot total rainfall from TTTs
per_ttt_plot=False      # plot percentage rainfall from TTTs (tot_ttt/tot_all)
rain_per_ttt_plot=False  # plot average rain per TTT day (rain composite)
comp_anom_ttt_plot=True  # plot rain per TTT as anom from long term daily mean for each month
comp_anom_ag_plot=False # plot comp anom with agtest on composite
comp_anom_cnt_plot=False     # plot count of the number of days above or below average
perc_ag=70              # show if this % or more days agree

# options to plot all - heavy pr
hvthrs=['0','10','25','50']
all_wet_cnt=False       # plot number of days over hvthr - either total or per mon depending on 'monmean'
all_wet_sum=False       # plot rainfall from days over hvthr - either total or per mon depending on 'monmean'
aper_wet_cnt=False       # plot % of days that are over hvthr
aper_wet_sum=False       # plot % of precip which is falling on days over hvthr

# options to plot ttt - heavy pr
ttt_wet_cnt=False        # plot number of  days over hvthr - either total or per mon depending on 'monmean'
tper_wet_cnt=False       # % of days over hvthr contributed by TTTs
per_tttd_wet=False       # plot % of TTT days which have precip over this threshold
ttt_wet_sum=False       # plot rainfall from days over hvthr - either total or per mon depending on 'monmean'
tper_wet_sum=False       # % of precip from days over hvthr contributed by TTTs


under_dayof='dayof'     # if "dayof" plots all rain on TTT days
                        #   if "under" plots rain under TTTs (based on blobs)
monmean='day'           # to control the output - is there averaging?
                        # 'day' is daily mean
                        # 'mon' is monthly mean
                        # 'tot' is total
nTTTlab=False            # labels each plot with # or % of TTTs

freecol=False           # free colour bar
refkey='0'              # 0 or all
#doms=['All']
doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected


bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
prdir=bkdir+"precip_figs/"

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

        ### Open rain data
        globp = 'pr'
        if dset == 'noaa':
            raindset = 'trmm'
            rainmod = 'trmm_3b42v7'
            rmoddct = dsetdict.dset_deets[raindset][rainmod]
            runits = rmoddct['prtimeunit']
            rcal = rmoddct['calendar']
            if testfile:
                rys = rmoddct['testfileyr']
            else:
                rys = rmoddct['yrfname']
            if testyear:
                rbeginatyr = rmoddct['testyr']
            else:
                rbeginatyr = rmoddct['startyr']
        else:
            raindset = dset
            rainmod = name
            rmoddct = moddct
            runits = units
            rcal = cal
            rys = ys
            rbeginatyr = beginatyr

        rainname = rmoddct['prname']
        rainfile = bkdir + raindset + "/" + rainmod + "."+globp+".day.mean." + rys + ".nc"
        print rainfile

        rainout = mync.open_multi(rainfile, globp, rainmod, \
                                  dataset=raindset, subs=subrain)

        rdim = len(rainout)
        if rdim == 5:
            rain, rtime, rlat, rlon, rdtime = rainout
        elif rdim == 6:
            rain, rtime, rlat, rlon, rlev, rdtime = rainout
            rain = np.squeeze(rain)
        else:
            print 'Check number of levels in ncfile'
        rdtime[:, 3] = 0

        ### Select data to run
        ### If testfile run on all days available
        if testfile:
            rain = rain[:, :, :];
            rtime = rtime[:];
            rdtime = rdtime[:]
        else:
            ### Find starting timestep
            start = rmoddct['startdate']
            ystart = int(start[0:4]);
            mstart = int(start[5:7]);
            dstart = int(start[8:10])
            if rcal == "360_day":
                startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
                beginday = ((int(rbeginatyr)) * 360) + 1
                daysgap = beginday - startday + 1
            else:
                startd = date(ystart, mstart, dstart)
                begind = date(int(rbeginatyr), 01, 01)
                daysgap = (begind - startd).days
            rain = rain[daysgap:, :, :];
            rtime = rtime[daysgap:];
            rdtime = rdtime[daysgap:]
        if testyear:
            if rcal == "360_day":
                rain, rdtime, rtime = rain[:360, :, :], rdtime[:360], rtime[:360]
            else:
                rain, rdtime, rtime = rain[:365, :, :], rdtime[:365], rtime[:365]

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
            if len(doms) == 3:
                kw, ke = stats.spatialsubset(s,False,cutlon=40.) # events west and east of 40E
                keys=[ks,kw,ke]
            elif len(doms) == 4:
                k1, ktmp = stats.spatialsubset(s, False, cutlon=37.5)
                k2, k3 = stats.spatialsubset(s,ktmp,cutlon=67.5)
                keys = [ks, k1, k2, k3]

            ### Plot rainmaps
            prbase=prdir+dset+"/"
            my.mkdir_p(prbase)
            mapsuf = seasopt+'_'+subrain+'_'+dset+'_'+name+'_'+thre_str+'_key'+refkey+'_4'+monmean
            if testfile or testyear:
                testq=True
            else:
                testq=freecol


            # Plots for all days ...

            if allplot:
                if t==0:
                    print 'Plotting all rain'
                    msklist=ap.gridrainmap_season(s,ks,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='tot_all',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=mapsuf,\
                                                  savefig=True,test=testq)

            if all_cnt_anom:
                if t == 0:
                    print 'Plotting count of all days with positive anom'
                    msklist = ap.gridrainmap_season(s, ks, rain, rlat, rlon, rdtime, units, cal, season=seasopt, \
                                                    key=dset + '-olr-0-' + refkey, ptype='all_cnt', mmean=monmean, \
                                                    under_of=under_dayof, figdir=prbase, file_suffix=mapsuf, \
                                                    savefig=True, test=testq)

            # Loop thresholds
            for h in range(len(hvthrs)):
                hstr=hvthrs[h]
                hvthr=float(hstr)
                print 'Plotting for heavy threshold '+hstr
                hvsuf=mapsuf+'_hvthr'+hstr


                if all_wet_cnt:
                    if t == 0:
                        print 'Plotting count of days over hvthr'
                        msklist = ap.gridrainmap_season(s, ks, rain, rlat, rlon, rdtime, units, cal, season=seasopt, \
                                                        key=dset + '-olr-0-' + refkey, ptype='all_wet_cnt', mmean=monmean, \
                                                        under_of=under_dayof, figdir=prbase, file_suffix=hvsuf, \
                                                        savefig=True,test=testq,labels=nTTTlab, heavy=hvthr)

                if all_wet_sum:
                    if t == 0:
                        print 'Plotting sum of pr from days over hvthr'
                        msklist = ap.gridrainmap_season(s, ks, rain, rlat, rlon, rdtime, units, cal, season=seasopt, \
                                                        key=dset + '-olr-0-' + refkey, ptype='all_wet_sum', mmean=monmean, \
                                                        under_of=under_dayof, figdir=prbase, file_suffix=hvsuf, \
                                                        savefig=True,test=testq,labels=nTTTlab,heavy=hvthr)

                if aper_wet_cnt:
                    if t == 0:
                        print 'Plotting % of days that are over hvthr'
                        msklist = ap.gridrainmap_season(s, ks, rain, rlat, rlon, rdtime, units, cal, season=seasopt, \
                                                        key=dset + '-olr-0-' + refkey, ptype='aper_wet_cnt', mmean=monmean, \
                                                        under_of=under_dayof, figdir=prbase, file_suffix=hvsuf, \
                                                        savefig=True, test=testq, labels=nTTTlab, heavy=hvthr)
                if aper_wet_sum:
                    if t == 0:
                        print 'Plotting % rainfall that falls on days over hvthr'
                        msklist = ap.gridrainmap_season(s, ks, rain, rlat, rlon, rdtime, units, cal, season=seasopt, \
                                                        key=dset + '-olr-0-' + refkey, ptype='aper_wet_sum', mmean=monmean, \
                                                        under_of=under_dayof, figdir=prbase, file_suffix=hvsuf, \
                                                        savefig=True, test=testq, labels=nTTTlab, heavy=hvthr)


            # Loop domains
            for do in range(len(doms)):
                domname=doms[do]
                print "Running on domain: "+domname
                eventkeys=keys[do]

                newsuf=mapsuf+'_'+domname

                if tot_ttt_plot:
                    print 'Plotting all rain from TTTs'
                    msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='tot_ttt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab)

                if per_ttt_plot:
                    print 'Plotting percentage rain from TTTs'
                    msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='per_ttt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab)

                if rain_per_ttt_plot:
                    print 'Plotting ave rain per TTT day'
                    msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='rain_per_ttt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq)

                if comp_anom_ttt_plot:
                    print 'Plotting composite rainfall anomalies'
                    msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ttt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab)

                if comp_anom_ag_plot:
                    newnewsuf=newsuf+'_agthr'+str(perc_ag)
                    print 'Plotting composite rainfall anomalies with ag test'
                    msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ag',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=newnewsuf,\
                                                  savefig=True,test=testq,labels=nTTTlab, agthresh=perc_ag)


                if comp_anom_cnt_plot:
                    msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_cnt',mmean=monmean,\
                                                  under_of=under_dayof,figdir=prbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq)

                # Loop thresholds
                for h in range(len(hvthrs)):
                    hstr = hvthrs[h]
                    hvthr = float(hstr)
                    print 'Plotting for TTT rain for heavy threshold ' + hstr
                    newhvsuf = mapsuf + '_hvthr' + hstr+'_' + domname


                    if ttt_wet_cnt:
                        print 'Plotting number of days over hvthr contributed by TTTs'
                        msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                      key=dset+'-olr-0-'+refkey,ptype='ttt_wet_cnt',mmean=monmean,\
                                                      under_of=under_dayof,figdir=prbase,file_suffix=newhvsuf,\
                                                      savefig=True,test=testq, heavy=hvthr)

                    if tper_wet_cnt:
                        print 'Plotting % of days over hvthr contributed by TTTs'
                        msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                      key=dset+'-olr-0-'+refkey,ptype='tper_wet_cnt',mmean=monmean,\
                                                      under_of=under_dayof,figdir=prbase,file_suffix=newhvsuf,\
                                                      savefig=True,test=testq, heavy=hvthr)

                    if ttt_wet_sum:
                        print 'Plotting hvthr day precip contributed by TTTs'
                        msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                      key=dset+'-olr-0-'+refkey,ptype='ttt_wet_sum',mmean=monmean,\
                                                      under_of=under_dayof,figdir=prbase,file_suffix=newhvsuf,\
                                                      savefig=True,test=testq, heavy=hvthr)

                    if tper_wet_sum:
                        print 'Plotting % hvthr day precip contributed by TTTs'
                        msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                      key=dset+'-olr-0-'+refkey,ptype='tper_wet_sum',mmean=monmean,\
                                                      under_of=under_dayof,figdir=prbase,file_suffix=newhvsuf,\
                                                      savefig=True,test=testq, labels=nTTTlab, heavy=hvthr)

                    if per_tttd_wet:
                        print 'Plotting % of TTTs with precip over hvthr'
                        msklist=ap.gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cal,season=seasopt,\
                                                      key=dset+'-olr-0-'+refkey,ptype='per_tttd_wet',mmean=monmean,\
                                                      under_of=under_dayof,figdir=prbase,file_suffix=newhvsuf,\
                                                      savefig=True,test=testq, labels=nTTTlab, heavy=hvthr)
