# Plotting wrapper
# to plot
# ....omega associated with TTTs
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
# but at present time have focused on noaa so the model stuff might not work well


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
import MetBot.dimensions_dict as dim_exdict
import MetBot.find_saddle as fs
import mpl_toolkits.basemap as bm

### Running options
sub="SA"
subvar="SA_TRMM"
#subvar="SA_CONT"
#subvar="UM_FOC"
seasopt="coreseason"    # options: coreseason, dryseason, fullseason
testyear=False           # To use output from a test
testfile=False           # Uses a test file with short period
                        # (testfile designed to be used together with testyear
                        # ..but testyear can be used seperately)
threshtest=False         # Option to run on thresholds + and - 5Wm2 as a test

# options to plot ttt...
comp_anom_ag_plot=True   # plot comp anom with agtest on composite
comp_anom_cnt_plot=True     # plot count of the number of days above or below average
perc_ag=80              # show if this % or more days agree

under_dayof='dayof'     # if "dayof" plots all rain on TTT days
                        #   if "under" plots rain under TTTs (based on blobs)

freecol=False           # free colour bar
refkey='0'              # 0 or all
#doms=['All']
doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected
varlist=['omega']
v=0 # for now don't loop variables, just one
levsel=True
if levsel:
    choosel=['500'] # can add a list
else:
    choosel=['1']
l=0 # for now don't loop variables, just one


bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
vardir=bkdir+varlist[v]+"_figs/"

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

    if dset != 'cmip5':
        levc = int(choosel[l])
    else:
        levc = int(choosel[l]) * 100

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

        # Get details - OLR
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

        ### Location for input & outputs
        indir=bkdir+dset+"/"
        outdir=indir+name+"/"
        if testyear: outdir=outdir+'test/'
        else: outdir=outdir
        my.mkdir_p(outdir)
        outsuf=outdir+name+'_'

        ### Open variable data
        globv = varlist[v]
        if dset == 'noaa' and globv != 'olr':
            if globv == 'pr':
                ds4noaa = 'trmm'
                mod4noaa = 'trmm_3b42v7'
            else:
                ds4noaa = 'ncep'
                mod4noaa = 'ncep2'
            dset2 = ds4noaa
            name2 = mod4noaa
        else:
            dset2 = dset
            name2 = name

        gmoddct=dsetdict.dset_deets[dset2][name2]
        gvname=gmoddct[globv+'name']
        if testfile:
            gys=gmoddct['testfileyr']
        else:
            if globv != 'omega' and globv != 'q' and globv != 'gpth':
                gys = gmoddct['yrfname']
            else:
                gys = gmoddct['fullrun']
        if testyear:
            gbeginatyr=gmoddct['testyr']
        else:
            gbeginatyr = gmoddct['startyr']
        gunits = gmoddct[globv+'timeunit']
        gcal = gmoddct['calendar']
        dimdict = dim_exdict.dim_deets[globv][dset2]
        glatname = dimdict[1]
        glonname = dimdict[2]

        # Open raw file
        rawfile = bkdir + dset2 + '/' + name2 + \
                  '.' + globv + '.day.mean.' + gys + '.nc'

        print 'Opening file: '+rawfile

        if levsel:
            ncout = mync.open_multi(rawfile, globv, name2, \
                                    dataset=dset2, subs=sub, levsel=levc)
        else:
            ncout = mync.open_multi(rawfile, globv, name2, \
                                    dataset=dset2, subs=sub)
        ndim = len(ncout)
        if ndim == 5:
            rawdata, gtime, glat, glon, gdtime = ncout
        elif ndim == 6:
            rawdata, gtime, glat, glon, glev, gdtime = ncout
            rawdata = np.squeeze(rawdata)
        else:
            print 'Check number of dims in ncfile'
        gdtime[:, 3] = 0

        ### Select data to run
        ### If testfile run on all days available
        if testfile:
            rawdata = rawdata[:, :, :];
            gtime = gtime[:];
            gdtime = gdtime[:]
        else:
            ### Find starting timestep
            start = gmoddct['startdate']
            ystart = int(start[0:4]);
            mstart = int(start[5:7]);
            dstart = int(start[8:10])
            if gcal == "360_day":
                startday = (ystart * 360) + ((mstart - 1) * 30) + dstart
                beginday = ((int(gbeginatyr)) * 360) + 1
                daysgap = beginday - startday + 1
            else:
                startd = date(ystart, mstart, dstart)
                begind = date(int(gbeginatyr), 01, 01)
                daysgap = (begind - startd).days
            rawdata = rawdata[daysgap:, :, :];
            gtime = gtime[daysgap:];
            gdtime = gdtime[daysgap:]
        if testyear:
            if gcal == "360_day":
                rawdata, gdtime, gtime = rawdata[:360, :, :], gdtime[:360], gtime[:360]
            else:
                rawdata, gdtime, gtime = rawdata[:365, :, :], gdtime[:365], gtime[:365]

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

            ### Plot varmaps
            gbase=vardir+dset+"/"
            my.mkdir_p(gbase)
            mapsuf = seasopt+'_'+subvar+'_'+dset+'_'+name+'_'+dset2+'_'+name2+'_'+thre_str+'_key'+refkey
            if testfile or testyear:
                testq=True
            else:
                testq=freecol


            # Loop domains
            for do in range(len(doms)):
                domname=doms[do]
                print "Running on domain: "+domname
                eventkeys=keys[do]

                newsuf=mapsuf+'_'+domname

                if comp_anom_ag_plot:
                    newnewsuf=newsuf+'_agthr'+str(perc_ag)
                    print 'Plotting composite anomalies with ag test'
                    msklist=ap.gridvarmap_season(s,eventkeys,globv,rawdata,glat,glon,gdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ag',\
                                                  under_of=under_dayof,figdir=gbase,file_suffix=newnewsuf,\
                                                  savefig=True,test=testq, agthresh=perc_ag)


                if comp_anom_cnt_plot:
                    msklist=ap.gridvarmap_season(s,eventkeys,globv,rawdata,glat,glon,gdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_cnt',\
                                                  under_of=under_dayof,figdir=gbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq)
