# Initial playing with SOMs

import sys,os
cwd=os.getcwd()
sys.path.append(cwd+'/..')
import MetBot.SynopticAnatomy as sy
import MetBot.MetBlobs as blb
import MetBot.dset_dict as dsetdict
import MetBot.find_saddle as fs


### Running options
testyear=True           # To use output from a test
testfile=True           # Uses a test file with short period
threshtest=True         # Option to run on thresholds + and - 5Wm2 as a test

bkdir=cwd+"/../../../CTdata/metbot_multi_dset/"

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

            ### Location for input & outputs
            indir=bkdir+dset+"/"
            outdir=indir+name+"/"
            if testyear: outdir=outdir+'test/'
            else: outdir=outdir
            outsuf = outdir + name + '_'

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

            ###  Open synop file
            syfile=outsuf+thre_str+'_'+dset+'-OLR.synop'
            s = sy.SynopticEvents((),[syfile],COL=False)