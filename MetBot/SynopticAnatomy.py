''' SynoptcAnatomy_fixing.py: A module of MetBot package\
Rainfall is an dictionary attribute of each event'''
import MetBlobs as blb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import mytools as my
import datetime
import cPickle
from time import time as timer
from time import sleep as sleep
import os
import shutil
import glob
figdir='/home/neil/work/computervision/metblobs/synopfigs/'

# CLASSES
class SynopticEvents:
    '''SynopticEvents(mbskeylist,mbslist,sbst='SA',hrwindow=49,COL=False)

    This class can contains synoptic events with associated meteorological
    fields.
    Needs a reference metblob set then stores associations to other metblob sets
    with criteria of a time window and space window
    Three initialization options exist:
    SynopticEvents(metblobskeylist,[metblobs1,...,metblobsN])
    where each metblobsN passed is tuple = (mbs, mbt, ch)

    SynopticEvents(metblobskeylist,[metblobs1,...,metblobsN])
    where files were created my Metblobs.mbsave: they open tuples as above

    SynopticEvents([],PreviouslySavedSEfile.synop)
    '''
    def __init__(self,mbskeylist,mbslist,sbst='SA',hrwindow=49,COL=False):
        fstr="hrtime, Label, Angle, cX, cY, minX, maxX, minY, maxY, area, cb.m10, cb.m01, cb.m11, cb.m20, cb.m02, cb.u11, cb.u20, cb.u02, cb.n11, cb.n20, cb.n02, cb.p1, cb.p2, Circularity"
        self.fields=[code.strip() for code in fstr.split(',')]

        domains={'SA': ((-60.0,0.0),(0.0,80.0)),
                 'Fauch08': ((-40.0,-15.0),(7.5,70.0),),
                'SH': ((-90.0,0.0),(0,360)),}
        vlt, vln = domains[sbst][0], domains[sbst][1]
        self.vaxis = vln[0],vln[1],vlt[0],vlt[1]

        self.mbskeys = tuple(mbskeylist)
        self.hrwindow = hrwindow
        self.blobs={}
        self.tracks={}
        self.trackshr={}
        self.trackscX={}
        self.trackscY={}
        self.events={}

        refm = mbslist[0]
        if isinstance(refm,tuple):
            for i in xrange(len(mbskeylist)):
                self.blobs[mbskeylist[i]] = {'mbs': mbslist[i][0],
                                             'mbt': mbslist[i][1],
                                             'ch' : mbslist[i][2]}
        elif isinstance(refm,str):
            #f, extens = refm.split('.')
            extens = refm.split('.')[-1]
            if extens == 'mbs':
                for i in xrange(len(mbskeylist)):
                    mbs, mbt, ch = blb.mbopen(mbslist[i])
                    self.blobs[mbskeylist[i]] = {'mbs': mbs,
                                                 'mbt': mbt,
                                                 'ch' : ch}
            elif extens == 'synop':
                tmo=timer()
                self.mbskeys, self.hrwindow, self.flagkey, self.blobs,\
                self.tracks, self.trackshr, self.trackscX, self.trackscY,\
                self.events = self.open(refm)
                print "Opening time:",(timer()-tmo),"s"
                if len(self.events)>0:
                    ks = self.events.keys()
                    print "Computing trkarrs of",len(ks),"events...";
                    if not COL:addtrkarrs(self)
                    elif COL:  addCOLs(self)
                    u = self.uniqueevents()

                print "Completion time:",(timer()-tmo),"s"
        self.flagkey = self.mbskeys[0] # only here because for synop call.. 
        # ...mbskeylist[0] won't exist at start

    def save(self,fname):
        print "Saving Synoptic Events to", fname
        pickf = open(fname,'w')
        dumptuple=(self.mbskeys, self.hrwindow,self.flagkey, self.blobs,\
          self.tracks, self.trackshr, self.trackscX, self.trackscY, self.events)
        cPickle.dump(dumptuple,pickf)
        pickf.close()

    def open(self,fname):
        '''flagkey, blobs, tracks, events = open(fname)'''
        print "Opening", fname
        pickf = open(fname,'r')
        mbskeys, hrwindow, flagkey, blobs, tracks, trackshr, trackscX,\
                                         trackscY, events = cPickle.load(pickf)
        pickf.close()
        return mbskeys, hrwindow, flagkey, blobs, tracks, trackshr, trackscX, trackscY, events

    # HELPER METHODS
    def ifld(self,codestr):
        '''Function for getting the index of metblobs dataset properties'''
        return self.fields.index(codestr)

    def __spheredist__(self,lon1, lat1, lon2, lat2):
        '''Calculated great-circle distances in metres'''
        R=6367442.76
        # Determine proper longitudinal shift
        dlon=np.abs(lon2-lon1)
        dlon=np.where(dlon>=180,360-dlon,dlon)
        #  Convert Decimal degrees to radians.
        deg2rad=np.pi/180
        lat1=lat1*deg2rad
        lat2=lat2*deg2rad
        dlon=dlon*deg2rad
        #
        #  Compute the distances
        t1 = np.sin(dlon) * np.cos(lat2)
        t2 = np.sin(lat2) * np.cos(lat1)
        t3 = np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
        dist=R*np.arcsin(np.sqrt(t1**2 + (t2-t3)**2))

        return dist

    def __rainamounts__(self,event,raindata,rainkey,heavy=20.):
        '''Calculates rainfall statistics for under cloud band footprints
        properties: rmn,rmx,wetness,heavyness,outrmn,outrmx,outwetness, rsum, hsum'''
        heavythresh=heavy
        rain,date,xypts = raindata
        erain=[]
        # STATION RAINFALL
        if rain.ndim==2:
            #print 'Station dataset ',rainkey
            for t in xrange(len(event.trkdtimes)):
                ix = my.ixdtimes(date,[event.trkdtimes[t,0]],\
                              [event.trkdtimes[t,1]],[event.trkdtimes[t,2]],[0])
                if len(ix)>1: print 'We have a problem'
                elif len(ix)==0:
                    if t==0:print 'No time match in',rainkey,\
                                  'for event starting on',event.trkdtimes[t]
                    continue
                chmask = points_inside_poly(xypts,\
                         event.blobs[event.refkey]['ch'][event.trk[t]])
                wetmask = rain[:,ix] >= 1.0
                rmask = chmask & wetmask.squeeze()
                heavymask = rain[:,ix] >= heavythresh
                hmask = chmask & heavymask.squeeze()
                r=np.ma.MaskedArray(rain[:,ix],mask=~rmask)
                h=np.ma.MaskedArray(rain[:,ix],mask=~hmask)
                if np.any(rmask):
                    # RAIN FALL INSIDE OLR COUNTOUR
                    rmn=r.mean()
                    rmx=r.max()
                    rsum=r.sum()
                    hsum=h.sum()
                    wetness=len(np.where(rmask)[0])/\
                            float(len(np.where(chmask)[0]))
                    heavyness=len(np.where(hmask)[0])/\
                              float(len(np.where(rmask)[0]))
                    # RAIN FALL OUTSIDE OLR COUNTOUR
                    ormask= ~chmask & wetmask.squeeze()
                    outr=np.ma.MaskedArray(rain[:,ix],mask=~ormask)
                    ormn=outr.mean()
                    ormx=outr.max()
                    owetness=len(np.where(ormask)[0])/\
                             float(len(np.where(~chmask)[0]))
                else:
                    rmn=np.NaN;rmx=np.NaN;wetness=np.NaN;heavyness=np.NaN
                    ormn=np.NaN;ormx=np.NaN;owetness=np.NaN
                erain.append((rmn,rmx,wetness,heavyness,ormn,ormx,owetness,rsum,hsum))
        # GRIDDED RAINFALL
        elif rain.ndim==3:
            #print 'Gridded rainfall dataset:',rainkey
            for t in xrange(len(event.trkdtimes)):
                ix = my.ixdtimes(date,[event.trkdtimes[t,0]],\
                              [event.trkdtimes[t,1]],[event.trkdtimes[t,2]],[0])
                #ix = my.ixdtimes(date,[event.trkdtimes[t,0]],\
                #              [event.trkdtimes[t,1]],[event.trkdtimes[t,2]],[event.trkdtimes[t,3]])
                if len(ix)>1: print 'We have a problem'
                elif len(ix)==0:
                    if t==0: print 'No time match in',rainkey,event.trkdtimes[t]
                    continue
                chmask = my.poly2mask(xypts[0],xypts[1],\
                                  event.blobs[event.refkey]['ch'][event.trk[t]])
                wetmask = rain[ix,:,:] >= 1.0
                rmask = chmask & wetmask
                heavymask = rain[ix,:,:] >= heavythresh
                hmask = chmask & heavymask
                r=np.ma.MaskedArray(rain[ix,:,:],mask=~rmask)
                h=np.ma.MaskedArray(rain[ix,:,:],mask=~hmask)
                if np.any(rmask):
                    # RAIN FALL INSIDE OLR COUNTOUR
                    rmn=r.mean()
                    rmx=r.max()
                    rsum=r.sum()
                    hsum=np.nansum(h)
                    wetness=len(np.where(rmask.ravel())[0])/\
                            float(len(np.where(chmask.ravel())[0]))
                    heavyness=len(np.where(hmask.ravel())[0])/\
                              float(len(np.where(rmask.ravel())[0]))
                    # RAIN FALL OUTSIDE OLR COUNTOUR
                    ormask= ~chmask & wetmask
                    outr=np.ma.MaskedArray(rain[ix,:,:],mask=~ormask)
                    ormn=outr.mean()
                    ormx=outr.max()
                    owetness=len(np.where(ormask.ravel())[0])/\
                             float(len(np.where(~chmask.ravel())[0]))
                else:
                    rmn=np.NaN;rmx=np.NaN;wetness=np.NaN;heavyness=np.NaN
                    ormn=np.NaN;ormx=np.NaN;owetness=np.NaN;rsum=np.NaN;hsum=np.NaN
                erain.append((rmn,rmx,wetness,heavyness,ormn,ormx,owetness,rsum,hsum))

        event.rainfall[rainkey] = np.asarray(erain)

    # BLOB MATCHING AND TRACKING

    def __mbsmatch__(self,mbs,ch,maxdist=1000e3):
        '''Matches metblobs[time] in blobs[now] to blobs[now+1]
        Matching very arbitrary: Uses criteria of only matching centroid within
                                 maxdist. There is however some commented
                                 code within that could allow dist criterion 
                                 to be 1.5xAveRadius of previous blob
        Returns: Array with one entry for each blob as such
        [date, label, nextdate, nextlabel, nextIndex, previousIndex]
        where array= -1  no matches'''

        # Get averaged radius of blob
        blobradii=[]
        icX, icY = self.ifld('cX'), self.ifld('cY')

        print "Calculating average blob convex hull radii..."
        for i in xrange(len(ch)):
            chlon, chlat = ch[i][:,0], ch[i][:,1]
            cx, cy = mbs[i, icX:icY+1]
            cx = np.tile(cx,(len(chlon),)) ; cy = np.tile(cy,(len(chlon),))
            #avedist = np.mean(self.__spheredist__(cx, cy, chlon, chlat))
            #if avedist > maxdist: avedist=maxdist
            avedist=maxdist
            ### arbitrary to use crit. as within 1.5xAveRadius of previous blob
            blobradii.append(avedist*1.0) 

        ### Match blobs by blobradii distance criterion
        mt=np.unique(mbs[:,0])
        ### This gives the time interval between timesteps
        dt = mt[1:]-mt[:-1];dt = dt.min() 
        matchnext=np.zeros((mbs.shape[0],6))
        matchnext[:]=-1
        matchnext[:,0:2]=mbs[:,0:2]
        print "Matching blobs..."
        ### Loop through each timestep
        for i in xrange(len(mt)-1):
            #print "Matching blobs found on: %d, (tstep=%d/%d)"\
            #        %(int(mt[i]), i, len(mt))
            ### Indices of blobs at timestep
            inow = np.where(mbs[:,0] == mt[i])[0]
            ### and following timestep
            inxt = np.where(mbs[:,0] == mt[i+1])[0]
            ### This makes sure don't match tracks beyond consecutive timesteps
            ### so gives an entry which only contains one blob. This implies
            ### a track will be built with only one entry
            if (mt[i+1]-mt[i]) > dt:
                for inw in inow:
                    matchnext[inw,:] = mbs[inw,0], mbs[inw,1], -1, -1, -999,\
                                       matchnext[inw,5]
                continue
            # For each blob centroid at now timestep
            for inw in inow:
                ### compare distance to all blob centroids at next timestep
                cx, cy = mbs[inw, icX:icY+1]
                cx = np.tile(cx,(len(inxt),)) ; cy = np.tile(cy,(len(inxt),))
                xnxt, ynxt = mbs[inxt,icX], mbs[inxt,icY]
                d2nxt = self.__spheredist__(cx, cy, xnxt, ynxt)
                centroid_moves_west = (xnxt - cx) > 0
                centroid_close_enough = d2nxt < blobradii[inw]
                m = np.where(centroid_close_enough & centroid_moves_west)[0]

                if len(m) > 0:
                    ### Get closest match if many pass criterion
                    ixbest = d2nxt.argmin()
                    ixmatch = inxt[ixbest]
                    matchnext[inw,:] = mbs[inw,0], mbs[inw,1], mbs[ixmatch,0],\
                                       mbs[ixmatch,1], ixmatch, matchnext[inw,5]
                    matchnext[ixmatch,5] = inw
                elif len(m)==0:
                    matchnext[inw,:] = mbs[inw,0], mbs[inw,1], -1, -1, -999,\
                                       matchnext[inw,5]
                    ### This is a key and important change!!! from what was
                    #matchnext[inw,:] = mbs[inw,0], mbs[inw,1], -1, -1, -1,\
                    #                   matchnext[inw,5]
                      
        matches = np.int32(matchnext)

        return matches

    def __blobtracks__(self,mbs,ch,trklen=2,maxdist=1000e3):
        '''Match subsequent blobs and create Tracks data structure
        Works because matches structure for each blob is:
        [date, label, nxtdate, nxtlabel, nxtIndex, prevIndex]

        Returns dictionary of tracks'''
        matches = self.__mbsmatch__(mbs,ch,maxdist=maxdist)
        # a view into just the Index fields of matches
        print 'Building tracks index dictionary...'
        m = matches[:,4:]
        tracks={}
        # Since track start occurs with no prevIndex but has a nxtIndex
        nxtNonxt = (m[:,0]>0) | (m[:,0]==-999)
        istart=np.where((m[:,1]==-1) & nxtNonxt )[0]
        for ix in istart:
            trk=[ix,]
            inxt = m[ix,0]
            if inxt == -999:
                dummy="this track will only have one position"
            else:
                while inxt>0:
                    trk.append(inxt)
                    inxt=m[inxt,0]
            if len(trk)>=trklen:
                trklabel = matches[ix,0]*100+matches[ix,1]
                tracks[trklabel] = trk

        return tracks

    def __trkmbsfield__(self,fieldstr,tracks,blobs):
        '''Return given fields from blobs for tracks'''
        ixf = self.ifld(fieldstr)
        field = tracks.copy()
        for k in tracks.keys():
            field[k]=list(blobs[tracks[k],ixf])
        return field

    def __tracks2blobs__(self,target,tracks,tracksmbs):
        '''Associated tracks with candidate/target blobs drawn from track's
           base blobs.
        Premise: Find which blobs exist in which tracks and associate
        Returns: array(iblob,blobkey,trackkey)'''

        print "Associating flag blobs with candidate tracks..."
        ihrt,ilab = self.ifld('hrtime'),self.ifld('Label')
        icx, icy = self.ifld('cX'), self.ifld('cY')
        startfilter = 30*24*100  # the time in label is *100
        # not sure if this is compatible with other datasets

        refhr, reflab = target[:,ihrt], target[:,ilab]
        refcx, refcy = target[:,icx], target[:,icy]
        refkey = np.int32(refhr)*100 + reflab
        ks=tracks.keys();ks.sort()
        ks=np.asarray(ks)
        trkshr = self.__trkmbsfield__('hrtime',tracks,tracksmbs)
        trkscX = self.__trkmbsfield__('cX',tracks,tracksmbs)
        trkscY = self.__trkmbsfield__('cY',tracks,tracksmbs)

        targettuple=[]
        for i in xrange(len(refhr)):
            rhr, rcX, rcY, rk = refhr[i], refcx[i], refcy[i], refkey[i]
            #print rhr, rlab
            ### First reduce size of search to only tracks 
            ### which start +-30days from target day
            ifltr = np.where(np.abs(ks-rk) < startfilter)[0]
            if len(ifltr) < 2: print len(ifltr)
            for j in ks[ifltr]:
                try:
                    idx = trkshr[j].index(rhr)
                    ### Now test if these blob features share a centroid: 
                    ### to avoid minor precision errors, we just check they
                    ### are very close, not exact.
                    diffCx = np.abs(trkscX[j][idx]-rcX)
                    diffCy = np.abs(trkscY[j][idx]-rcY)
                    if diffCx<1.0 and diffCy<1.0:
                        targettuple.append((i,rk,j))
                        #print 'yes this worked for',j
                    else:
                        dummy='dummy'
                        #print j,'FAILED'
                except ValueError:
                    continue

        return np.int32(np.asarray(targettuple))


    # PLOTTTING METHODS
    def __blobpatch__(self,xverts,yverts):
        verts=np.hstack((xverts[:,np.newaxis],yverts[:,np.newaxis]))
        #verts = verts-verts.mean(0)
        Path = mpath.Path
        moves = ','.join(['Path.LINETO' for i in xrange(verts.shape[0]-2)])
        exec('codes = [Path.MOVETO,'+moves+', Path.CLOSEPOLY]')
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path,facecolor='b',alpha=0.7)
        #plt.text(verts.mean(0)[0],verts.mean(0)[1] , "Cloudband",\
        #          ha="center", size=14)
        ax = plt.axes()
        ax.add_patch(patch)


    def __plotblob__(self,blobentry,chentry,patch=False,bmap=False,**kwargs):
        '''Simple blob plot
        Can pass basemapobj if desire, otherwise standard plt.plot performed'''
        if not bmap: bmap = plt
        plotprops=[]
        for k in kwargs:
            plotprops.append(k+"='"+kwargs[k]+"'")
        if len(plotprops) == 0: plotprops = ["'b--'","linewidth=2"]
        plotprops = ','.join(plotprops)

        cX, cY = blobentry[self.ifld('cX')], blobentry[self.ifld('cY')]
        xverts = np.append(chentry[:,0],chentry[0,0])
        yverts = np.append(chentry[:,1],chentry[0,1])

        bmap.plot(cX,cY,'kx',markersize=12)
        exec('bmap.plot(xverts,yverts,'+plotprops+')')
        if patch:
            self.__blobpatch__(xverts,yverts)



    def __plottrack__(self, trkmbskey,trkkey, bmap=False,**kwargs):
        '''Simple track plot'''
        trk = self.tracks[trkmbskey][trkkey]
        coast=True
        if not bmap:
            coast=False
            bmap = plt
        plotprops=[]
        for k in kwargs:
            plotprops.append(k+"='"+kwargs[k]+"'")
        if len(plotprops) == 0: plotprops = ["'r-'","linewidth=2"]
        plotprops = ','.join(plotprops)

        trkcXs = self.blobs[trkmbskey]['mbs'][trk,self.ifld('cX')]
        trkcYs = self.blobs[trkmbskey]['mbs'][trk,self.ifld('cY')]
        # Plot blobs being tracked
        cnt=1
        for i in trk:
            if coast: bmap.drawcoastlines();bmap.drawcountries()
            self.__plotblob__(self.blobs[trkmbskey]['mbs'][i],\
                              self.blobs[trkmbskey]['ch'][i], bmap=bmap)
            exec('bmap.plot(trkcXs[:cnt],trkcYs[:cnt],'+plotprops+')')
            plt.axis(self.vaxis)
            tm=self.blobs[trkmbskey]['mbt'][i]
            humandate=' Date: %d-%02d-%02d %02dh00' %(tm[0],tm[1],tm[2],tm[3])
            plt.title(trkmbskey+" Trk: "+str(trkkey)+humandate)
            plt.draw()
            cnt+=1
            sleep(0.1)
            plt.clf()

    def seetracks(self,trkmbskey,bmap=False):
        '''Loop through and plot tracks for given mbskey'''
        #plt.figure();plt.show()
        trks = self.tracks[trkmbskey]
        keyslist = trks.keys();keyslist.sort()
        for k in keyslist:
            self.__plottrack__(trkmbskey,k, bmap=bmap)
            plt.clf()

    # MAIN ROUTINES OF CLASS SYNOPTIC EVENTS
    def buildtracks(self,keys=False,trklen=1,maxdist=1000e3):
        '''Will build tracks for all metblobs unless keys are specified'''
        self.tracks={}
        self.trackshr={}
        self.trackscX={}
        self.trackscY={}
        
        keylist=self.mbskeys[:]
        if keys:keylist=keys
        for k in keylist:
            print 'GETTING TRACKS FOR:', k
            # GET THE TRACKS
            self.tracks[k] = self.__blobtracks__(self.blobs[k]['mbs'],\
                             self.blobs[k]['ch'],trklen=trklen,maxdist=maxdist)
            # GET KEY PROPERTIES NEEDED IN buildevents below
            self.trackshr[k] = self.__trkmbsfield__('hrtime',self.tracks[k],\
                                                    self.blobs[k]['mbs'])
            self.trackscX[k] = self.__trkmbsfield__('cX',self.tracks[k],\
                                                    self.blobs[k]['mbs'])
            self.trackscY[k] = self.__trkmbsfield__('cY',self.tracks[k],\
                                                     self.blobs[k]['mbs'])
            print 'SUCCESS!\n'

    def buildevents(self,basetrkkey='noaa-olr-0-0',trklen=1,maxdist=1000e3):
        '''Builds list of Events
        Builds events using a prescribed base track specified by basetrkkey
        Default basetrkkey = "noaa-olr-0-all"
        but reason to believe better to use noaa-olr-0-0 which would require
        trklen in buildtracks to be set = 1
        '''
        self.events = {}
        if len(self.tracks)==0:
            print '''First need to build tracks...here we go...
            You should have provided trklen and maxdist keyword arguments to to
            clear of your result that follows...
            '''
            self.buildtracks(trklen=trklen,maxdist=maxdist)
        # MATCH REFERENCE BLOBTRACKS TO CANDIDATE DAYS
        ix_base = self.mbskeys.index(basetrkkey)
        flagmbs = self.blobs[self.flagkey]['mbs']
        candtracks = self.tracks[self.mbskeys[ix_base]]
        candblobs = self.blobs[self.mbskeys[ix_base]]['mbs']
        
        if len(candtracks)==0:
            print "Cannot proceed as no candidate base tracks available!"
            return False
        asstrks = self.__tracks2blobs__(flagmbs,candtracks, candblobs)
        self.asstrks=asstrks
        ixflagblob, flagtrackkeys = asstrks[:,0], asstrks[:,2]
        Bequeath = self
        count=1
        print 'Building events...'

        tms=timer()
        for k in np.unique(flagtrackkeys):
            tmk=timer()
            iflag = np.where(flagtrackkeys == k)[0]
            an_event = Event(Bequeath, k, ixflagblob[iflag],\
                             basetrkkey=basetrkkey,maxdist_otherblobs=3000e3)
            self.events[k]=an_event
            print 'Built event: %d/%d in %4.2f s'\
                   %(count,len(np.unique(flagtrackkeys)),(timer()-tmk))
            count += 1
        print "...that took",(timer()-tms)/60,"mins"
        uniques = self.uniqueevents()
        addtrkarrs(self)

    def addeventrain(self,rainkeys,type='station',heavy=20., \
                     datadir='/home/neil/sogehome/data/'):
        '''rainkeys can be 'wrc', 'trmm' but must be list
        type is either 'station' or 'grid'
        Would make more sense to put this partly in buildevents and the
        self.rainamounts function into a method of Event but for various 
        reason, not least, need to open and close rain datasets,
        it is done this way'''
        # STATION DATA
        ekeys = self.events.keys();ekeys.sort()
        if type=='station':
            for rainkey in rainkeys:
                print 'Adding rain from ',rainkey,'station data set'
                exec('rain, hrtime, lat, lon, dates = my.open%sstations(\
                      datadir=\"%s\")'\
                       %(rainkey.upper(),datadir))
                rainstations=(rain,dates,np.hstack((lon,lat)))
                for k in ekeys:
                    evnt = self.events[k]
                    self.__rainamounts__(evnt,rainstations,rainkey,heavy=heavy)
                del rain, dates, lat, lon, rainstations
        if type=='grid':
            for rainkey in rainkeys:
                print 'Adding rain from ',rainkey,'gridded data set'
                flist=glob.glob('%s/%s/%s.*.daily.SA.nc'\
                                 %(datadir,rainkey,rainkey) )
                flist.sort()
                rain, time, lat, lon, lev, dtime = \
                            my.OpenMultipleNC(flist,rainkey,'ncep2',sub='SA')
                dtime[:,3]=0
                raingrid=(rain,dtime,(lon,lat))
                for k in ekeys:
                    evnt=self.events[k]
                    self.__rainamounts__(evnt,raingrid,rainkey+'grid',\
                                         heavy=heavy)
                del rain, dtime, lat, lon, raingrid

    def addeventrain_any(self, raindata, lat, lon, dates,\
                         rainkeys, type='grid', heavy=20.):
        '''rainkeys can be 'wrc', 'trmm' but must be list
        type is either 'station' or 'grid'
        Would make more sense to put this partly in buildevents and the
        self.rainamounts function into a method of Event but for various
        reason, not least, need to open and close rain datasets,
        it is done this way'''
        # Edited to allow raindata to be set in wrapper
        ekeys = self.events.keys();
        ekeys.sort()
        if type == 'grid':
            for rainkey in rainkeys:
                print 'Adding rain from ', rainkey, 'gridded data set'
                dates[:, 3] = 0 # to set hrtime to 0
                raingrid = (raindata, dates, (lon, lat))
                for k in ekeys:
                    evnt = self.events[k]
                    self.__rainamounts__(evnt, raingrid, rainkey + 'grid', \
                                         heavy=heavy)
                del dates, lat, lon, raingrid

    def uniqueevents(self):
        '''Loop through all events and ensure none give same event but with
        different starting blobs'''
        print "Checking for uniqueness of events.."
        tm=timer
        ks = self.events.keys();ks.sort()
        uniques=[]
        for i in xrange(len(ks)):
            samevents=[ks[i]]
            e = self.events[ks[i]]
            trkarr = e.trkarrs[e.refkey]
            nxt=i+10
            if nxt>(len(ks)-1): nxt = len(ks)-1

            icnt=1
            while (i+icnt) < nxt:
                ecomp = self.events[ks[i+icnt]]
                otrkarr = ecomp.trkarrs[e.refkey]
                nmatches=0
                for tm, ib in trkarr[:,:2]:
                    if np.any((tm==otrkarr[:,0]) & (ib==otrkarr[:,1]))\
                    and ib !=-1:
                        nmatches += 1
                if nmatches>0: samevents.append((ks[i+icnt],nmatches))
                icnt+=1

            uniques.append(samevents)

        self.uniques = uniques
        return uniques

class Event(SynopticEvents):
    '''Object that contains all details of single synoptic events
       Might be clumsy, but inherits from SynopticEvents, so has Synoptic Events
       methods.
    '''

    def __init__(self,Inherits,reftrackkey,ixflags,maxdist_otherblobs=1000e3,\
                 sbst='SA',basetrkkey='noaa-olr-0-all'):
        #print 'Initiating event: ',reftrackkey
        self.fields = Inherits.fields
        self.blobs = Inherits.blobs
        self.tracks = Inherits.tracks
        self.trackshr = Inherits.trackshr
        self.trackscX = Inherits.trackscX
        self.trackscY = Inherits.trackscY
        self.mbskeys = Inherits.mbskeys
        self.mxdist = maxdist_otherblobs
        self.vaxis = Inherits.vaxis
         
        self.refkey = basetrkkey
        self.trkkey = reftrackkey
        self.trk = self.tracks[self.refkey][self.trkkey]
        self.ixflags = ixflags
        self.flagtimes = self.blobs[Inherits.flagkey]['mbs']\
                                   [ixflags,self.ifld('hrtime')]
        self.trktimes = self.blobs[self.refkey]['mbs']\
                        [self.trk,self.ifld('hrtime')]
        self.trkdtimes = self.blobs[self.refkey]['mbt'][self.trk,:]
        self.trkcX = self.blobs[self.refkey]['mbs'][self.trk,self.ifld('cX')]
        self.trkcY = self.blobs[self.refkey]['mbs'][self.trk,self.ifld('cY')]
        self.assoctrks = {}
        if len(self.mbskeys) > 2:
            self.ombskeys = self.mbskeys[2:]
            for k in self.ombskeys:
                self.assoctrks[k] = self.__othertracks__(self.trkkey,k,\
                                                         maxdist=self.mxdist)
        # Build eventarrays needed for further plotting and analysis
        self.trkarrs,refkey = self.__eventarray__()
        #self.trkarrs = {} # this also can get added later by calling self.__eventarray__
        self.rainfall = {} # this gets added in later outside of this class by SynopticEvents.addeventrain()


    def __othertracks__(self,flagtrackkey,k,maxdist=1000e3):
        '''Associate event track with other tracks
        Returns: array(trackkey, ntimeoverlap)'''
        ihrt, ilab = self.ifld('hrtime'), self.ifld('Label')
        icX, icY = self.ifld('cX'), self.ifld('cY')
        startfilter = 30*24*100  # the time in label is *100
        #.....not sure if the above is compatible with other datasets?
        # Get necessary properties of elements of flag track
        flagtrk = self.trk
        flaghrs = self.trktimes
        flagcX  = self.trkcX
        flagcY  = self.trkcY

        trkshr = self.trackshr[k]
        trkscX = self.trackscX[k]
        trkscY = self.trackscY[k]
        # Get array of keys of the other tracks and filter for +- 90 days
        # either side
        ks = trkshr.keys();ks.sort()
        ks=np.array(ks)
        # Get indices of tracks to search for associations
        ifltr = np.where(np.abs(ks-flagtrackkey) < startfilter)[0]
        if len(ifltr) < 2: print len(ifltr)

        assocs=[]
        for k in ks[ifltr]:
            ohrs, ocX, ocY = np.asarray(trkshr[k]),\
                             np.asarray(trkscX[k]), np.asarray(trkscY[k])
            cnt=0
            for i in xrange(len(flaghrs)):
                rhr, rcX, rcY = flaghrs[i], flagcX[i], flagcY[i]
                sphered = self.__spheredist__(rcX, rcY, ocX.mean(),ocY.mean())
                overlap = (ohrs[0] <= rhr) & (rhr <= ohrs[-1])
                if (sphered < maxdist) and overlap:
                    cnt += 1

            if cnt > 0: assocs.append((k,cnt))

        return assocs

    def __eventarray__(self):
        '''Method of Events to create an array that describes event
        array has dimensions time x properties x variable
        properties: hrtime, trkix, cX, cY
        variable: these are those contained in mbskeys

        Rainfall is attribute later to this object by simply gettting
        stations/gridpoints within the OLR cloud contour its a dict with keys
        for each rainfall dataset each key accesses an array time x values
        values are meanrain,maxrain,wetness,heavyness,outsidemean,
        outsidemax,outsidewetness
        '''
        # Get hour time for event from appearance of first metblob in any track
        mbskeyvals=self.refkey.split('-')
        self.ombskeys = self.mbskeys[2:]
        allhrtime=np.asarray(self.trktimes)
        for mbsk in self.ombskeys:
            for k, o in self.assoctrks[mbsk]:
                allhrtime=np.append(allhrtime,self.trackshr[mbsk][k])
        allhrtime=np.unique(allhrtime)
        shp = (len(allhrtime),4)
        trkarr=np.zeros(shp,dtype=np.float32);trkarr[:]=-1
        # BUILD TRACK ARRAY FOR REFERENCE/FLAG TRACKS: 
        # currently 'noaa-olr-0-all' which has tres==24
        allhrtm = allhrtime
        #if mbskeyvals[0]=='noaa':
            # round the time (in hours) to the nearest day?
            # not sure this does anything because the value seems to be zero?
            # for other datasets it might break things
            # ....because than allhrtm does not equal t below
        #    allhrtm = allhrtm - (allhrtm % 24)
        trkarr[:,0]=allhrtm
        cnt=0
        for t in self.trktimes:
            #print self.trktimes, allhrtime, allhrtm
            #print np.where(allhrtm==t)
            ix = np.where(allhrtm==t)[0][0]
            mbt = self.blobs[self.refkey]['mbt'][self.trk[cnt]]
            trkarr[ix:ix+4,:] = t,self.trk[cnt],self.trkcX[cnt],self.trkcY[cnt]
            cnt+=1

        ##### Build dictionary of tracks
        trkarrs={self.refkey:np.int32(trkarr[:,:-1])}
        trkarrs={self.refkey:trkarr}
        trkarrstime={self.refkey:trkarr[:,0]}

        # BUILD TRACK ARRAYS FOR OTHER ARRAYS
        for mbsk in self.ombskeys:
            # Get associated tracks for given mbskey
            #dset,vrb,lvl,descr=mbsk.split('-')
            otrks = self.assoctrks[mbsk]
            if len(otrks)==0: # Since requires existence of associated track
                print 'No associated',mbsk,'tracks for event',\
                       self.trkkey,'moving on...'
                trkarrs[mbsk], trkarrstime[mbsk] = {}, {}
                continue

            shp=(len(allhrtime),4,len(otrks))
            otrkarr=np.zeros(shp,dtype=np.float32);otrkarr[:]=-1
            #if mbsk=='COL': otrkarr[:,0,:]=allhrtm
            for i in xrange(len(otrks)): otrkarr[:,0,i]=allhrtime
            # Get sorted list of associated track keys
            trkks=[]
            for k,o in otrks: trkks.append(k)
            trkks.sort()

            trkcnt=0
            for k in trkks:
                otrk = self.tracks[mbsk][k]
                otrkhr = self.trackshr[mbsk][k]
                otrkcX = self.trackscX[mbsk][k]
                otrkcY = self.trackscY[mbsk][k]
                cnt=0
                for t in otrkhr:
                    ix = np.where(allhrtime==t)[0]
                    mbt = self.blobs[mbsk]['mbt'][otrk[cnt]]
                    otrkarr[ix,:,trkcnt] = t,otrk[cnt],otrkcX[cnt],otrkcY[cnt]
                    if mbsk=='COL':
                        ix = np.where(allhrtm==t)[0][0]
                        otrkarr[ix:ix+4,:,trkcnt] = t,otrk[cnt],\
                                                    otrkcX[cnt],otrkcY[cnt]
                    cnt+=1
                trkcnt+=1

            trkarrs[mbsk] = otrkarr
            if otrkarr.ndim==2:
                trkarrstime[mbsk]=otrkarr[:,0]
            if otrkarr.ndim==3:
                trkarrstime[mbsk]=otrkarr[:,0,0]
        return trkarrs, trkarrstime

# FUNCTIONS: EXTRA EVENT ATTRIBUTES
def gridmasks(event,lon,lat,trkkey=False,compress=True):
    '''This simple script produces grid masks for each contour associated with
        each event
    Usage: masks = event.gridmasks(lat,lon)
          lat and lon are the coordinates of the associated gridded data
          trkkey=False implies get the masks for all mbskeys
          trkkey='noaa-olr-0-all' implies compute the masks for only
                 that dataset
          compress=True means multi-track variables are going to be flattened
                        into one mask

    Returns: masks (dict) mask['mbs-keys-you-specified'] = {}
    '''
    e = event
    if trkkey==False:
        mbskeys = e.trkarrs.keys()
        mbskeys.sort()
    elif isinstance(trkkey,str):
        mbskeys = [trkkey]
    elif isinstance(trkkey,list):
        mbskeys = trkkey
    # MAIN LOOP OF METHOD
    masks = {}
    for mbk in mbskeys:
        if isinstance(e.trkarrs[mbk],dict):
            print 'No track array for',mbk,'in event',e.trkkey
            continue
        trkarr = np.int32(e.trkarrs[mbk])
        if trkarr.ndim==2:
            ntimes = trkarr.shape[0]
            maskarr = np.bool8(np.ones((ntimes,len(lat),len(lon))))
            for t in xrange(ntimes):
                if trkarr[t,1] != -1:
                    maskarr[t,:,:] = my.poly2mask(lon,lat,event.blobs[mbk]['ch'][trkarr[t,1]],invert=True)

        elif trkarr.ndim==3:
            ntimes, ntrks = trkarr.shape[0],trkarr.shape[2]
            maskarr = np.bool8(np.ones((ntimes,len(lat),len(lon),ntrks)))
            for t in xrange(ntimes):
                for n in xrange(ntrks):
                    if trkarr[t,1,n] != -1:
                        maskarr[t,:,:,n] = my.poly2mask(lon,lat,e.blobs[mbk]['ch'][trkarr[t,1,n]],invert=True)
            if compress:
                maskarr_c = maskarr[:,:,:,0]
                for z in xrange(ntrks):
                    maskarr_c = maskarr_c & maskarr[:,:,:,z]
                maskarr = maskarr_c

        masks[mbk] = maskarr
    return masks

def coordsmasks(event,xy_points,trkkey=False):
    '''This simple script produces masks for coordinates within each contour
       associated with each event
    NOTE: This produces masks for station RAINFALL DATA
    Usage: masks = event.gridmasks(event,xy_points)
          xy_points are the coordinates of the associated list of points
          xy_points[:,0] is lon, xy_points[:,1] is lat
          trkkey=False implies get the masks for all mbskeys
          trkkey='noaa-olr-0-all' implies compute the masks for only
                 that dataset

    Returns: masks (dict) mask['mbs-keys-you-specified'] = {}
    '''
    e = event
    if trkkey==False:
        mbskeys = e.trkarrs.keys()
    elif isinstance(trkkey,str):
        mbskeys = [trkkey]
    elif isinstance(trkkey,list):
        mbskeys = trkkey
    # MAIN LOOP OF METHOD
    masks = {}
    npts = xy_points.shape[0]
    for mbk in mbskeys:
        trkarr = np.int32(e.trkarrs[mbk])
        if trkarr.ndim==2:
            ntimes = trkarr.shape[0]
            maskarr = np.bool8(np.zeros((ntimes,npts)))
            for t in xrange(ntimes):
                if trkarr[t,1] != -1:
                    maskarr[t,:] = points_inside_poly(xy_points,\
                                               e.blobs[mbk]['ch'][trkarr[t,1]])

        elif trkarr.ndim==3:
            ntimes, ntrks = trkarr.shape[0],trkarr.shape[2]
            maskarr = np.bool8(np.zeros((ntimes,n,ntrks)))
            for t in xrange(ntimes):
                for n in xrange(ntrks):
                    if trkarr[t,1,n] != -1:
                        maskarr[t,:,n] = points_inside_poly(xy_points,\
                                               e.blobs[mbk]['ch'][trkarr[t,1]])
        masks[mbk] = maskarr
    return masks

def addtrkarrs(SEobj):
    '''This is generally performed at file load time in 
       SynopticAnatomy.SynopticEvents.__init__() however can be called outside
       of this.
       Unecessary to call this if already plan on calling 
       SynopticAnatomy.addCOLs()
    '''
    ks = SEobj.events.keys()
    print "Computing trkarrs of",len(ks),"events..."
    for k in ks:
        refkey = SEobj.events[k].refkey
        SEobj.events[k].trkarrs, SEobj.events[k].trkarrstime =\
                                               SEobj.events[k].__eventarray__()


def addCOLs(SEobj,mintrklen=2,fname='/home/neil/phd/cutofflowtracks_alicefavre/ncep2tracks.txt'):
    '''Adds cut-off low tracks from Alice Favre's tracks'''
    COLs, COLhrtime, COLdtime = my.readincutoffs(fname=fname)
    nl = COLs.shape[0];dum=np.tile(0,(nl,1))
    SEobj.blobs['COL'] =\
    {'mbs': np.hstack((COLhrtime[:,np.newaxis],COLs[:,1][:,np.newaxis],\
                       COLs[:,2][:,np.newaxis],COLs[:,4][:,np.newaxis],\
                       COLs[:,3][:,np.newaxis])),\
    'mbt': COLdtime}
    trkkeys = np.unique(COLs[:,0])
    SEobj.tracks['COL']={}
    for k in trkkeys:
        keypart='%02d' %(k)
        trk=list(np.int32(np.where(COLs[:,0]==k)[0]))
        # CONSTRUCT HRTIME+FEATURE_NO KEY
        key=int(str(int(COLhrtime[trk[0]]))+keypart[-2:])
        if len(trk)>=mintrklen:
            SEobj.tracks['COL'][key] = trk
    newk='COL'
    SEobj.trackshr[newk] = SEobj.__trkmbsfield__('hrtime',SEobj.tracks[newk],\
                                                 SEobj.blobs[newk]['mbs'])
    SEobj.trackscX[newk] = SEobj.__trkmbsfield__('cX',SEobj.tracks[newk],\
                                                 SEobj.blobs[newk]['mbs'])
    SEobj.trackscY[newk] = SEobj.__trkmbsfield__('cY',SEobj.tracks[newk],\
                                                 SEobj.blobs[newk]['mbs'])

    oldmbskeys=list(SEobj.mbskeys);oldmbskeys.append('COL')
    SEobj.mbskeys=tuple(oldmbskeys)

    for k in SEobj.events.keys():
        e=SEobj.events[k];e.mbskeys=SEobj.mbskeys
        e.assoctrks[newk] = e.__othertracks__(e.trkkey,newk,maxdist=1000e3)
        e.trkarrs, e.trkarrstime = e.__eventarray__()

def addRWB(SEobj,mintrklen=1,wb='AWB',pv=1.5,isk=340,fname='/home/neil/work/rwb/fromThando/'):
    '''Adds cut-off low tracks from Thando Ndarana's RWB events'''
    rwb, rwbhrtime, rwbdtime = my.readin_rwb(wb=wb,pv=pv,isk=isk,fdir=fname)
    nl = rwb.shape[0]
    rwbkeys = np.arange(1,nl+1,1)
    SEobj.blobs[wb] =\
    {'mbs': np.hstack((rwbhrtime[:,np.newaxis],rwbkeys[:,np.newaxis],\
                       rwb[:,2][:,np.newaxis],rwb[:,3][:,np.newaxis],\
                       rwb[:,4][:,np.newaxis])),\
     'mbt': rwbdtime}
    trkkeys = rwbkeys
    SEobj.tracks[wb]={}
    for k in trkkeys:
        keypart='%04d' %(k)
        #trk=list(np.int32(np.where(rwb[:,1]==k)[0]))
        trk=[k-1]
        print trk
        # CONSTRUCT HRTIME+FEATURE_NO KEY
        key=int(str(int(rwbhrtime[trk[0]]))+keypart)
        print key
        if len(trk)>=mintrklen:
            SEobj.tracks[wb][key] = trk
    newk=wb
    SEobj.trackshr[newk] = SEobj.__trkmbsfield__('hrtime',SEobj.tracks[newk],\
                                                  SEobj.blobs[newk]['mbs'])
    SEobj.trackscX[newk] = SEobj.__trkmbsfield__('cX',SEobj.tracks[newk],\
                                                  SEobj.blobs[newk]['mbs'])
    SEobj.trackscY[newk] = SEobj.__trkmbsfield__('cY',SEobj.tracks[newk],\
                                                 SEobj.blobs[newk]['mbs'])

    oldmbskeys=list(SEobj.mbskeys);oldmbskeys.append(wb)
    SEobj.mbskeys=tuple(oldmbskeys)

def points_inside_poly(points,poly):
    thepoly=Path(poly)
    boolean=thepoly.contains_points(points)
    return boolean

    #for k in SEobj.events.keys():
        #e=SEobj.events[k];e.mbskeys=SEobj.mbskeys
        #e.assoctrks[newk] = e.__othertracks__(e.trkkey,newk,maxdist=1000e3)
        #e.trkarrs, e.trkarrstime = e.__eventarray__()
