# Plotting wrapper
# to plot
# ....qflux associated with TTTs
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
import scipy.interpolate as spi
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

threshtest=False         # Option to run on thresholds + and - 5Wm2 as a test

# options to plot ttt...
comp_anom_ttt_plot=True  # plot ave OLR on TTT days as anom from long term daily mean for each month
comp_anom_ag_plot=True   # plot comp anom with agtest on composite
perc_ag=70              # show if this % or more days agree

under_dayof='dayof'     # if "dayof" plots all rain on TTT days
                        #   if "under" plots rain under TTTs (based on blobs)

freecol=False           # free colour bar
refkey='0'              # 0 or all
swapd='ncep'
#doms=['All']
doms=['All','nCont','nMada','nOcea'] # doms for TTT days selected
varlist=['qflux']
v=0 # for now don't loop variables, just one
levsel=True
if levsel:
    choosel=['850'] # can add a list
else:
    choosel=['1']
l=0 # for now don't loop variables, just one
interp=True
int_res=3.5

bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
vardir=bkdir+varlist[v]+"_figs/"
my.mkdir_p(vardir)

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
        ys=moddct['yrfname']
        beginatyr = moddct['startyr']
        units = moddct['olrtimeunit']
        cal = moddct['calendar']

        ### Location for input & outputs
        indir=bkdir+dset+"/"
        outdir=indir+name+"/"
        my.mkdir_p(outdir)
        outsuf=outdir+name+'_'

        ### Open variable data
        variable = varlist[v]
        globv=variable
        if dset == 'noaa' and globv != 'olr':
            if globv == 'pr':
                ds4noaa = 'trmm'
                mod4noaa = 'trmm_3b42v7'
            else:
                if swapd=='era':
                    ds4noaa='era'
                    mod4noaa='erai'
                elif swapd=='ncep':
                    ds4noaa = 'ncep'
                    mod4noaa = 'ncep2'
            dset2 = ds4noaa
            name2 = mod4noaa
        else:
            dset2 = dset
            name2 = name

        print "Running on " + variable
        if variable == 'wind':
            globv1 = 'u'
            globv2 = 'v'
        elif variable == 'qflux':
            globv1 = 'u'
            globv2 = 'v'
            globv3 = 'q'

        gmoddct=dsetdict.dset_deets[dset2][name2]
        if globv == 'qflux':
            gys = gmoddct['fullrun']
        else:
            print 'Need to specify years'
        gbeginatyr = gmoddct['startyr']
        gunits = gmoddct[globv1+'timeunit']
        gcal = gmoddct['calendar']
        dimdict = dim_exdict.dim_deets[globv1][dset2]
        glatname = dimdict[1]
        glonname = dimdict[2]

        vnamedict_u = globv1 + 'name'
        vnamedict_v = globv2 + 'name'
        if variable == 'qflux':
            vnamedict_q = globv3 + 'name'

        varstr_u = gmoddct[vnamedict_u]
        varstr_v = gmoddct[vnamedict_v]
        if variable == 'qflux':
            varstr_q = gmoddct[vnamedict_q]

        # Open raw file
        rawfile_u = bkdir  + dset2 + '/' + name2 + \
                    '.' + globv1 + '.day.mean.' + gys + '.nc'
        rawfile_v = bkdir + dset2 + '/' + name2 + \
                    '.' + globv2 + '.day.mean.' + gys + '.nc'

        if variable == 'qflux':
            rawfile_q = bkdir + dset2 + '/' + name2 + \
                        '.' + globv3 + '.day.mean.' + gys + '.nc'

        print 'Opening ' + rawfile_u
        print 'and corresponding file: ' + rawfile_v
        if variable == 'qflux':
            print 'and q file:' + rawfile_q

        if levsel:
            ncout_u = mync.open_multi(rawfile_u, globv1, name2, \
                                      dataset=dset2, subs=sub, levsel=levc)
            ncout_v = mync.open_multi(rawfile_v, globv2, name2, \
                                      dataset=dset2, subs=sub, levsel=levc)

            if variable == 'qflux':
                ncout_q = mync.open_multi(rawfile_q, globv3, name2, \
                                          dataset=dset2, subs=sub, levsel=levc)

        else:
            ncout_u = mync.open_multi(rawfile_u, globv1, name2, \
                                      dataset=dset2, subs=sub)
            ncout_v = mync.open_multi(rawfile_v, globv2, name2, \
                                      dataset=dset2, subs=sub)

            if variable == 'qflux':
                ncout_q = mync.open_multi(rawfile_q, globv3, name2, \
                                          dataset=dset2, subs=sub)

        ndim = len(ncout_u)

        if ndim == 5:

            rawdata_u, gtime, glat, glon, gdtime = ncout_u
            rawdata_v, gtime, glat, glon, gdtime = ncout_v

            if variable == 'qflux':
                rawdata_q, gtime, glat, glon, gdtime = ncout_q

        elif ndim == 6:

            rawdata_u, gtime, glat, glon, glev, gdtime = ncout_u
            rawdata_u = np.squeeze(rawdata_u)

            rawdata_v, gtime, glat, glon, glev, gdtime = ncout_v
            rawdata_v = np.squeeze(rawdata_v)

            if variable == 'qflux':
                rawdata_q, gtime, glat, glon, glev, gdtime = ncout_q
                rawdata_q = np.squeeze(rawdata_q)

        else:
            print 'Check number of dims in ncfile'

        gdtime[:, 3] = 0

        nlat = len(glat)
        nlon = len(glon)
        nsteps = len(rawdata_u[:, 0, 0])

        if interp:
            print "Interpolating data to a " + str(int_res) + " grid"
            minlon = glon[0]
            maxlon = glon[nlon - 1]
            minlat = glat[0]
            maxlat = glat[nlat - 1]

            newlon = np.arange(minlon, maxlon + 2.5, int_res)
            newlat = np.arange(maxlat, minlat + 2.5, int_res)
            # newlat=newlat[::-1]

            nlon = len(newlon)
            nlat = len(newlat)

            prodata_u = np.zeros((nsteps, nlat, nlon), dtype=np.float32)
            prodata_v = np.zeros((nsteps, nlat, nlon), dtype=np.float32)

            if variable == 'qflux':
                prodata_q = np.zeros((nsteps, nlat, nlon), dtype=np.float32)

            # Get rid of nans
            nonan_u = np.nan_to_num(rawdata_u)
            nonan_v = np.nan_to_num(rawdata_v)

            if variable == 'qflux':
                nonan_q = np.nan_to_num(rawdata_q)

            for step in range(nsteps):
                Interpolator_u = spi.interp2d(glon, glat, nonan_u[step, :, :], kind='linear')
                Interpolator_v = spi.interp2d(glon, glat, nonan_v[step, :, :], kind='linear')

                prodata_u[step, :, :] = Interpolator_u(newlon, newlat)
                prodata_v[step, :, :] = Interpolator_v(newlon, newlat)

                if variable == 'qflux':
                    Interpolator_q = spi.interp2d(glon, glat, nonan_q[step, :, :], kind='linear')
                    prodata_q[step, :, :] = Interpolator_q(newlon, newlat)

        else:

            prodata_u = rawdata_u
            prodata_v = rawdata_v

            if variable == 'qflux':
                prodata_q = rawdata_q

            newlon = lon
            newlat = lat

        # If qflux then multiply winds with humidity
        if variable == 'qflux':
            print "Multiplying winds by q..."

            qu = prodata_u * prodata_q
            qv = prodata_v * prodata_q

            prodata_u = qu
            prodata_v = qv

        ### Get thresholds and loop
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
            testq=freecol


            # Loop domains
            for do in range(len(doms)):
                domname=doms[do]
                print "Running on domain: "+domname
                eventkeys=keys[do]

                newsuf=mapsuf+'_'+domname

                if comp_anom_ttt_plot:
                    print 'Plotting composite anomalies'
                    msklist=ap.gridvectmap_season(s,eventkeys,globv,prodata_u,prodata_v,newlat,newlon,gdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ttt',\
                                                  under_of=under_dayof,figdir=gbase,file_suffix=newsuf,\
                                                  savefig=True,test=testq, agthresh=perc_ag)


                if comp_anom_ag_plot:
                    newnewsuf=newsuf+'_agthr'+str(perc_ag)
                    print 'Plotting composite anomalies with ag test'
                    msklist=ap.gridvectmap_season(s,eventkeys,globv,prodata_u,prodata_v,newlat,newlon,gdtime,units,cal,season=seasopt,\
                                                  key=dset+'-olr-0-'+refkey,ptype='comp_anom_ag',\
                                                  under_of=under_dayof,figdir=gbase,file_suffix=newnewsuf,\
                                                  savefig=True,test=testq, agthresh=perc_ag)