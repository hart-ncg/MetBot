# Wrapper to output a list of the timesteps of OLR file should be used for SOM analysis
# The timesteps are those with TTCBs (currently based on all days within events)

import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import numpy as np
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs
import MetBot.mynetcdf as mync
import MetBot.mytools as my


### Running options
testyear=True           # To use output from a test
testfile=True           # Uses a test file with short period
if testfile:
    print "Running on testfiles"
threshtest=True         # Option to run on thresholds + and - 5Wm2 as a test

### Get directory
bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"
sub="SA"

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
            mnames=['noaa']
        nmstr=str(nmod)

        for m in range(nmod):
            name=mnames[m]
            mcnt=str(m+1)
            print 'Running on ' + name
            print 'This is model '+mcnt+' of '+nmstr+' in list'

            # Get details
            moddct=dsetdict.dset_deets[dset][name]
            if testfile:
                ys = moddct['testfileyr']
            else:
                ys = moddct['yrfname']
            if testyear:
                beginatyr = moddct['testyr']
            else:
                beginatyr = moddct['startyr']

            ### Location for input & outputs
            indir=bkdir+dset+"/"
            outdir=indir+name+"/"
            infile = indir + name + ".olr.day.mean." + ys + ".nc"
            if testyear: outdir=outdir+'test/'
            else: outdir=outdir
            outsuf = outdir + name + '_'

            ### Open OLR data - select region but all dates are fine for now
            print "Opening " + infile
            v = dset + "-olr-0-0"
            daset, globv, lev, drv = v.split('-')
            ncout = mync.open_multi(infile, globv, name,\
                                    dataset=dset, subs=sub)
            ndim = len(ncout)
            if ndim == 5:
                olr, time, lat, lon, dtime = ncout
            elif ndim == 6:
                olr, time, lat, lon, lev, dtime = ncout
                olr = np.squeeze(olr)
            else:
                print 'Check number of dims in ncfile'

            ### Find dates of TTT events
            "Getting dates of TTT events"

            ### Get thresh
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

            ### Select events
            ks = s.events.keys()
            ks.sort()

            ### Count number of events
            count_all=str(int(len(ks)))
            print "Total CB events ="+str(count_all)

            ### Get a list of dates - as np.array (nstep,4) yyyy, mm, dd, hh
            ### I wanted to check the indices using two methods for now
            refkey = s.mbskeys[0]
            edts = []
            indices_m1=[]
            for k in ks:
                e = s.events[k]
                dts = s.blobs[refkey]['mbt'][e.ixflags]
                for dt in range(len(dts)):
                    edts.append(dts[dt])
                    # Check if it matches dtime
                    ix = my.ixdtimes(dtime, [dts[dt][0]], \
                                     [dts[dt][1]], [dts[dt][2]], [0])
                    indices_m1.append(ix)
            edts = np.asarray(edts)
            indices_m1 = np.asarray(indices_m1)

            ### Compare all dates with the list of TTT dates
            # 2nd method
            str_edts = [str(edts[x]) for x in range(len(edts[:, 0]))]
            str_dtime = [str(dtime[x]) for x in range(len(dtime[:, 0]))]

            indices_m2=[]
            for i in range(len(str_edts)):
                indices_m2.append(str_dtime.index(str_edts[i]))

            ### Check methods get the same result
            if len(indices_m1)==len(indices_m2):
                print "Dates of CBs detected for "+str(len(indices_m1))+" days"
            else:
                print "Error in detecting number of days with CBs"
            for c in range(len(edts)):
                if indices_m1[c]!=indices_m2[c]:
                    print "Error in detecting days with CBs"

            ### Sel OLR (time,lat,lon) for TTCB days
            olrsel=olr[indices_m1,:,:]

            ### Print list of indices
            outtxt="../../SOM_test/txtfile_tstep/"+dset+"."+name+".thresh_"+thnames[t]+".tstep_ttcb.txt"
            print "Outputing tsteps to txtfile "+outtxt
            txtfile = open(outtxt, "w")
            for i in range(len(indices_m1)):
                txtfile.write(","+str(indices_m1[i][0]))
            txtfile.close()