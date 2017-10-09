# AdvancedPlots.py
import SynopticAnatomy as sy
import numpy as np
import matplotlib.pyplot as plt
import EventStats as stats
#import SynopticClassify as sc
import SynopticPlot as syp
import RainStats as rs
import mytools as my
import MetBlobs as blb
from time import sleep as sleep

### SUBPLOT FUNCTIONS
def spatiofreq2_season(s,lat,lon,yrs,eventkeys,figno=1,season='coreseason',\
     key='noaa-olr-0-all',flagonly=False,file_suffix='test',savefig=False,\
     fontdict=False):
    '''spatiofreq2_season(s,lat,lon,yrs,figno=1,season='coreseason',\
                          flagonly=False,file_suffix='test')
    Produces subplots of cloud-band gridpoint count by month
    '''
    if not fontdict:fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbkl=key.split('-')
    if mbkl[0]=='noaa':dclim=(1,9,1)
    if mbkl[0]=='hadam3p':dclim=(1,11,1)
    if mbkl[0]=='um':dclim=(1,9,1)
    if mbkl[0]=='cmip5':dclim=(1,9,1)

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season=='dryseason':mns=[4,5,6,7,8,9]
    elif isinstance(season,list):
        mns=season
    m, f = blb.SAfrBasemap(lat[4:-7],lon[3:-2],drawstuff=True,prj='cyl',\
           fno=figno,rsltn='l',fontdict=fd)
    if len(mns)==12:plt.close();plt.figure(figsize=[12,10])
    elif len(mns)==6:plt.close();plt.figure(figsize=[13,8])
    cnt=1
    msklist=[]
    for mn in mns:
        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)
        if flagonly:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,\
                    clim=(1.,9.,1.),month=mn,flagonly=True,fontdict=fd)
        else:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,clim=dclim,\
                                      month=mn,flagonly=False,fontdict=fd)
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
        my.xtickfonts();my.ytickfonts()
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(file_suffix+'_flagonly.png',dpi=150)
        else: plt.savefig(file_suffix+'.png',dpi=150)

    return msklist

def spatiofreq3_season(s,lat,lon,yrs,eventkeys,figno=1,season='coreseason',\
     key='noaa-olr-0-all',res='noaa',dclim=(1,7,1),flagonly=False,file_suffix='test',savefig=False,\
     fontdict=False):
    '''spatiofreq3_season(s,lat,lon,yrs,figno=1,season='coreseason',\
                          flagonly=False,file_suffix='test')
    Produces subplots of cloud-band gridpoint count by month
    spatiofreq3 similar to spatiofreq2 but with option by RJ to fix res to noaa
    '''
    if not fontdict:fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbkl=key.split('-')

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season=='dryseason':mns=[4,5,6,7,8,9]
    elif isinstance(season,list):
        mns=season

    if mbkl[0]=='noaa':
        m, f = blb.SAfrBasemap(lat[6:-8],lon[3:-2],drawstuff=True,prj='cyl',fno=figno,rsltn='l')
    else:
        if res=='native':
            m, f = blb.SAfrBasemap(lat[13:-27],lon[6:-5],drawstuff=True,prj='cyl',fno=figno,rsltn='l')
        else:
            m, f = blb.SAfrBasemap(lat[6:-8], lon[3:-2], drawstuff=True, prj='cyl', fno=figno,
                                   rsltn='l')  # lat and lon specified by number - diff for UM at native res

    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])


    if len(mns)==12:plt.close();plt.figure(figsize=[12,10])
    elif len(mns)==6:plt.close();plt.figure(figsize=[13,8])
    cnt=1
    msklist=[]
    for mn in mns:
        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)
        if flagonly:
            allmask=stats.spatiofreq3(m,s,lat,lon,yrs,eventkeys,\
                    clim=dclim,month=mn,flagonly=True,fontdict=fd)
        else:
            allmask=stats.spatiofreq3(m,s,lat,lon,yrs,eventkeys,clim=dclim,\
                                      month=mn,flagonly=False,fontdict=fd)
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
        my.xtickfonts();my.ytickfonts()
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(file_suffix+'_res-'+res+'_flagonly_spatiofreq.png',dpi=150)
        else: plt.savefig(file_suffix+'_res-'+res+'_spatiofreq.png',dpi=150)

    return msklist

def spatiofreq2_seasonanoms(s,lat,lon,yrs,eventkeys,msklist,figno=1,\
    season='coreseason',key='noaa-olr-0-all',flagonly=False,\
    file_suffix='test',savefig=False,fontdict=False):
    '''spatiofreq2_season(s,lat,lon,yrs,figno=1,season='coreseason',\
                          flagonly=False,file_suffix='test')
    Produces subplots of cloud-band gridpoint count by month
    '''
    if not fontdict: fd = {'fontsize':14,'fontweight':'bold'}
    else: fd=fontdict
    mbkl=key.split('-')
    if mbkl[0]=='noaa':dclim=(-3.,3,.5)
    if mbkl[0]=='hadam3p':dclim=(-5.,5,1.)

    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season
    m, f = blb.SAfrBasemap(lat[4:-7],lon[3:-2],drawstuff=True,prj='cyl',\
                           fno=figno,rsltn='l',fontdict=fd)
    if len(mns)==12:plt.close();plt.figure(figsize=[12,10])
    elif len(mns)==6:plt.close();plt.figure(figsize=[13,8])
    cnt=1
    for mn in mns:
        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)
        if flagonly:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,clim=dclim,\
            meanmask=msklist[cnt-1],month=mn,flagonly=True,fontdict=fd)
        else:
            allmask=stats.spatiofreq2(m,s,lat,lon,yrs,eventkeys,fontdict=fd,\
            meanmask=msklist[cnt-1],clim=dclim,month=mn,flagonly=False)
        if cnt%2==0:
            syp.redrawmap(m,lns=True,resol='verylow',parallel=False,\
                          fontdict=fd,clr='0.5')
        else:
            syp.redrawmap(m,lns=True,resol='verylow',fontdict=fd,clr='0.5')
        cnt+=1
    plt.subplots_adjust(left=0.05,right=0.92,top=0.97,bottom=0.02,\
                        wspace=0.02,hspace=0.1)
    if savefig:
        if flagonly:
            plt.savefig(stats.figdir+file_suffix+'_cbcount_flagonly.png',\
                        dpi=150)
        else:
            plt.savefig(stats.figdir+file_suffix+'_cbcount.png',dpi=150)


def gridrainmap_season(s,eventkeys,rain,rlat,rlon,rdtime,units,cl,season='coreseason',key='noaa-olr-0-0',\
                       ptype='per_ttt',mmean='mon',under_of='dayof',figdir='test',file_suffix='test',\
                       savefig=False, test=True, labels=False,agthresh='perc_ag',heavy='0'):
    '''Produces subplots of ttt rainfall by month
    need to open the rain data with lon and lat and also the synop file
    see e.g. plot_ttt_precip_autothresh.py

    plot types
        tot_all - sum of total precip
        all_cnt - % of all days with +ve pr anoms

        all_wet_cnt - plot number of days over hvthr - either total or per mon depending on 'monmean'
        all_wet_sum - plot rainfall from days over hvthr - either total or per mon depending on 'monmean'
        aper_wet_cnt - % of days that are hvthr days
        aper_wet_sum - plot % of precip which is falling on hvthr days


        tot_ttt - sum of precip from TTTs
        per_ttt - percent of precip from TTTs
        rain_per_ttt - rain per ttt event
        comp_anom_ttt - rain per ttt as anom from monthly mean precip
        comp_anom_ag - as comp anom but with ag test
        comp_anom_cnt - % days with +ve or -ve anomalies

        ttt_wet_cnt- number of hvthr days - either total or per mon depending on 'monmean'
        tper_wet_cnt- % of hvthr days contributed by TTTs
        per_tttd_wet - % of TTT days which have precip over this threshold
        ttt_wet_sum- rainfall from hvthr days - either total or per mon depending on 'monmean'
        tper_wet_sum- % of hvthr day precip contributed by TTTs

    under_of -> "dayof" is rain on day of TTTs, "under" is rain under TTTs
    '''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    yrs = np.unique(rdtime[:,0])
    nys=len(yrs)

    # Get list of dates for these events
    # and if 'under' a list of chs
    edts = []
    if under_of == 'under':
        chs = []
    ecnt=1
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[key]['mbt'][e.ixflags]
        #dts=e.trkdtimes
        for dt in range(len(dts)):
            if ecnt==1:
                edts.append(dts[dt])
                if under_of=='under':
                    chs.append(e.blobs[key]['ch'][e.trk[dt]])
            else:
                tmpdt=np.asarray(edts)
                # Check if it exists already
                ix = my.ixdtimes(tmpdt, [dts[dt][0]], \
                     [dts[dt][1]], [dts[dt][2]], [0])
                if len(ix)==0:
                    edts.append(dts[dt])
                    if under_of == 'under':
                        chs.append(e.blobs[key]['ch'][e.trk[dt]])
            ecnt+=1
    edts = np.asarray(edts)
    edts[:, 3] = 0
    if under_of == 'under':
        chs=np.asarray(chs)
    print "Number of original TTT days found =  " + str(len(edts))

    # Get n lat and lon
    nlon=len(rlon)
    nlat=len(rlat)

    # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season == 'dryseason':mns = [4, 5, 6, 7, 8, 9]
    elif isinstance(season,list):
        mns=season
    print season
    print mns

    # Draw basemap
    m, f = blb.SAfrBasemap(rlat,rlon,drawstuff=True,prj='cyl',fno=1,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    cnt=1
    msklist=[]
    for mn in mns:
        print 'Plotting month '+str(mn)

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data (allmask)
        #All rain
        firstyear=yrs[0]
        lastyear=yrs[nys-1]
        raindat=np.where((rdtime[:,0]>=firstyear) & (rdtime[:,0]<=lastyear))
        rainperiod=rain[raindat,:,:]
        rainperiod=np.squeeze(rainperiod)
        newdates=rdtime[raindat]
        ix=np.where((newdates[:,1]==mn))
        rainmon=rainperiod[ix,:,:]
        datesmon=newdates[ix]
        ndays_mon=len(datesmon)
        rainmon=np.squeeze(rainmon)
        rainsum_all=np.nansum(rainmon,0)
        rainsum_monthly=rainsum_all/nys
        rainsum_daily=rainsum_all/len(datesmon)

        #Days over hvthr
        hvthr=float(heavy)
        wetdaycnt=np.zeros((nlat, nlon), dtype=np.float32)
        wetdaysum=np.zeros((nlat, nlon), dtype=np.float32)
        for i in range(nlat):
            for j in range(nlon):
                wet_ind = np.where(rainmon[:, i, j] > hvthr)[0]

                count = len(wet_ind)
                wet_sel = rainmon[wet_ind,i,j]
                wet_sum = np.nansum(wet_sel,0)

                wetdaycnt[i,j] = count
                wetdaysum[i,j] = wet_sum

        wetdays_p_mon=wetdaycnt/nys
        wetdays_mean=np.nanmean(wetdays_p_mon)
        wetdays_per=(wetdaycnt/ndays_mon)*100.0

        wetsum_p_mon=wetdaysum/nys
        wetsum_mean=np.nanmean(wetsum_p_mon)
        wetsum_per=(wetdaysum/rainsum_all)*100.0

        if ptype=='tot_all':
            if mmean=='mon':
                data4plot=rainsum_monthly
            elif mmean=='day':
                data4plot=rainsum_daily
            elif mmean=='tot':
                data4plot=rainsum_all

        elif ptype=='all_wet_cnt':
            if mmean=='tot':
                data4plot=wetdaycnt
            elif mmean=='mon':
                data4plot=wetdays_p_mon
            titstat=wetdays_mean

        elif ptype=='all_wet_sum':
            if mmean=='tot':
                data4plot=wetdaysum
            elif mmean=='mon':
                data4plot=wetsum_p_mon
            titstat=wetsum_mean

        elif ptype=='aper_wet_cnt':
            data4plot=wetdays_per

        elif ptype=='aper_wet_sum':
            data4plot=wetsum_per


        elif ptype=='all_cnt':

            all_anoms = np.zeros((ndays_mon, nlat, nlon), dtype=np.float32)
            for day in range(ndays_mon):
                this_anom = rainmon[day, :, :] - rainsum_daily
                all_anoms[day, :, :] = this_anom

            all_pos_pcent = np.zeros((nlat, nlon), dtype=np.float32)

            for i in range(nlat):
                for j in range(nlon):
                    count_p = len(np.where(all_anoms[:, i, j] > 0)[0])
                    perc_p = (float(count_p) / float(ndays_mon)) * 100
                    all_pos_pcent[i, j] = perc_p

            data4plot=all_pos_pcent

        else:

            #TTT rain
            ix2 = np.where((edts[:, 1] == mn))
            edatesmon = edts[ix2]
            if under_of=='under':
                chs_mon=chs[ix2]

            indices = []
            if under_of == 'under':
                chs_4rain = []
            for edt in range(len(edatesmon)):
                ix = my.ixdtimes(datesmon, [edatesmon[edt][0]], \
                                 [edatesmon[edt][1]], [edatesmon[edt][2]], [0])
                if len(ix) >= 1:
                    indices.append(ix)
                    if under_of=='under':
                        chs_4rain.append(chs_mon[edt])
            if len(indices) >= 2:
                indices = np.squeeze(np.asarray(indices))
                if under_of=='under':
                    chs_4rain=np.asarray(chs_4rain)
            else:
                indices = indices
                if under_of=='under':
                    chs_4rain=np.squeeze(np.asarray(chs_4rain))
            nttt_mon = len(indices)
            print "Number of TTT days found in rain dataset for mon " + str(mn) + " =  " + str(nttt_mon)

            if under_of == 'dayof':

                if nttt_mon == 0:
                    rainsum_ttt = np.zeros((nlat, nlon), dtype=np.float32)
                else:
                    rainsel = rainmon[indices, :, :]

                    if nttt_mon>=2:
                        rainsum_ttt=np.nansum(rainsel,0)
                    else:
                        rainsum_ttt=np.squeeze(rainsel)

                rainsum_ttt_monthly=rainsum_ttt/nys
                rainsum_ttt_daily=rainsum_ttt/ndays_mon

                rainperttt=rainsum_ttt/nttt_mon

                comp_anom = rainperttt - rainsum_daily

                # Wet days
                ttt_wetdaycnt = np.zeros((nlat, nlon), dtype=np.float32)
                ttt_wetdaysum = np.zeros((nlat, nlon), dtype=np.float32)
                for i in range(nlat):
                    for j in range(nlon):
                        wet_ind = np.where(rainsel[:, i, j] > hvthr)[0]
                        count = len(wet_ind)
                        wet_sel = rainsel[wet_ind, i, j]
                        wet_sum = np.nansum(wet_sel, 0)

                        ttt_wetdaycnt[i, j] = count
                        ttt_wetdaysum[i, j] = wet_sum

                ttt_wetdays_p_mon = ttt_wetdaycnt / nys
                ttt_wetdays_mean = np.nanmean(ttt_wetdays_p_mon)
                ttt_wetdays_per = (ttt_wetdaycnt / wetdaycnt)* 100.0
                per_ttt_wthis = (ttt_wetdaycnt / nttt_mon)  * 100.0

                ttt_wetsum_p_mon = ttt_wetdaysum / nys
                ttt_wetsum_mean = np.nanmean(ttt_wetsum_p_mon)
                ttt_wetsum_per = (ttt_wetdaysum / wetdaysum)* 100.0

                if ptype=='comp_anom_ag' or ptype=='comp_anom_cnt':
                    if nttt_mon >= 1:
                        anoms = np.zeros((nttt_mon, nlat, nlon), dtype=np.float32)
                        for day in range(nttt_mon):
                            this_anom = rainsel[day, :, :] - rainsum_daily
                            anoms[day, :, :] = this_anom

                        if ptype == 'comp_anom_ag':

                            anoms_signs = np.sign(anoms)
                            comp_signs = np.sign(comp_anom)


                            mask_zeros = np.zeros((nlat, nlon), dtype=np.float32)
                            for i in range(nlat):
                                for j in range(nlon):
                                    count = len(np.where(anoms_signs[:, i, j] == comp_signs[i, j])[0])
                                    perc = (float(count) / float(nttt_mon)) * 100
                                    if perc >= agthresh:
                                        mask_zeros[i, j] = 1
                                    else:
                                        mask_zeros[i, j] = 0


                        elif ptype=='comp_anom_cnt':

                            pos_pcent=np.zeros((nlat, nlon), dtype=np.float32)
                            neg_pcent=np.zeros((nlat, nlon), dtype=np.float32)
                            #zero_pcent=np.zeros((nlat, nlon), dtype=np.float32)

                            for i in range(nlat):
                                for j in range(nlon):
                                    count_p = len(np.where(anoms[:, i, j] > 0)[0])
                                    count_n = len(np.where(anoms[:, i, j] < 0)[0])
                                    #count_z = len(np.where(anoms[:, i, j] == 0)[0])


                                    perc_p = (float(count_p) / float(nttt_mon)) * 100
                                    perc_n = (float(count_n) / float(nttt_mon)) * 100
                                    #perc_z = (float(count_z) / float(nttt_mon)) * 100
                                    #print perc_z


                                    pos_pcent[i,j]=perc_p
                                    neg_pcent[i,j]=perc_n
                                    #zero_pcent[i,j]=perc_z

            elif under_of=='under':

                if nttt_mon == 0:
                    rainsum_ttt = np.zeros((nlat, nlon), dtype=np.float32)
                else:
                    rainsel = rainmon[indices, :, :]
                    ttt_rain_dates = datesmon[indices]
                    ndt=nttt_mon

                    masked_rain=np.ma.zeros((ndt,nlat,nlon),dtype=np.float32)
                    for rdt in range(ndt):
                        chmask = my.poly2mask(rlon,rlat,chs_4rain[rdt])
                        r=np.ma.MaskedArray(rainsel[rdt,:,:],mask=~chmask)
                        masked_rain[rdt,:,:]=r

                    if nttt_mon >= 2:
                        rainsum_ttt = np.ma.sum(masked_rain, 0) # maskd elements set to 0
                    else:
                        rainsum_ttt = np.ma.squeeze(masked_rain)

                rainsum_ttt_monthly = rainsum_ttt / nys
                rainsum_ttt_daily = rainsum_ttt / ndays_mon

                rainperttt = rainsum_ttt / nttt_mon


            if ptype=='tot_ttt':
                if mmean == 'mon':
                    data4plot = rainsum_ttt_monthly
                elif mmean == 'day':
                    data4plot = rainsum_ttt_daily
                elif mmean == 'tot':
                    data4plot = rainsum_ttt
            elif ptype=='per_ttt':
                if mmean == 'mon':
                    rainperc_ttt=(rainsum_ttt_monthly/rainsum_monthly)*100.0
                elif mmean == 'day':
                    rainperc_ttt = (rainsum_ttt_daily / rainsum_daily) * 100.0
                elif mmean == 'tot':
                    rainperc_ttt=(rainsum_ttt/rainsum_all)*100.0
                data4plot=rainperc_ttt
            elif ptype=='rain_per_ttt':
                data4plot=rainperttt
            elif ptype=='comp_anom_ttt':
                data4plot=comp_anom
            elif ptype=='comp_anom_ag':
                #data4plot=comp_anom[::3,::3]
                data4plot=comp_anom
            elif ptype=='comp_anom_cnt':
                data4plot=pos_pcent
            elif ptype=='ttt_wet_cnt':
                if mmean=='tot':
                    data4plot=ttt_wetdaycnt
                elif mmean=='mon':
                    data4plot=ttt_wetdays_p_mon
            elif ptype=='tper_wet_cnt':
                data4plot=ttt_wetdays_per
            elif ptype=='ttt_wet_sum':
                if mmean=='tot':
                    data4plot=ttt_wetdaysum
                elif mmean=='mon':
                    data4plot=ttt_wetsum_p_mon
            elif ptype=='tper_wet_sum':
                data4plot=ttt_wetsum_per
            elif ptype=='per_tttd_wet':
                data4plot=per_ttt_wthis

        newlon=rlon
        newlat=rlat

        #Plot
        plon,plat = np.meshgrid(newlon,newlat)

        # all plot cbars
        if ptype=='tot_all':
            if mmean=='tot':
                clevs=[0,800,1600,2400,3200,4000,4800,5600]
            elif mmean=='mon':
                clevs=[0,50,100,150,200,250,300,350]
            elif mmean=='day':
                clevs=[0.0,1.5,3.0,4.5,6.0,7.5,9.0,10.5,12.0]
            cm = plt.cm.YlGnBu
            cbar_lab='mm'
        elif ptype=='all_cnt':
            clevs=np.arange(0,50,5)
            cm = plt.cm.magma_r
            cbar_lab='%'

        # all hvthr plot cbars
        elif ptype=='all_wet_cnt' or ptype=='ttt_wet_cnt':
            cm = plt.cm.YlGnBu
            cbar_lab = 'days'
            if hvthr==0:
                if mmean == 'tot':
                    clevs=np.arange(0,400,40)
                elif mmean =='mon':
                    clevs=np.arange(0,22,2)
            elif hvthr==10:
                if mmean == 'tot':
                    clevs=np.arange(0,225,25)
                elif mmean =='mon':
                    clevs=np.arange(0,16,1)
            elif hvthr==25:
                if mmean == 'tot':
                    clevs=np.arange(0,60,5)
                elif mmean =='mon':
                    clevs=np.arange(0,5,0.5)
            elif hvthr==50:
                if mmean == 'tot':
                    clevs=np.arange(0,16,1)
                elif mmean =='mon':
                    clevs=np.arange(0,2,0.25)

        elif ptype=='all_wet_sum' or ptype=='ttt_wet_sum':
            cm = plt.cm.YlGnBu
            cbar_lab = 'mm'
            if hvthr==0:
                if mmean == 'tot':
                    clevs = [0, 800, 1600, 2400, 3200, 4000, 4800, 5600]
                elif mmean =='mon':
                    clevs = [0, 50, 100, 150, 200, 250, 300, 350]
            elif hvthr==10:
                if mmean == 'tot':
                    clevs=np.arange(0,3000,500)
                elif mmean =='mon':
                    clevs=np.arange(0,300,30)
            elif hvthr==25:
                if mmean == 'tot':
                    clevs=np.arange(0,2600,300)
                elif mmean =='mon':
                    clevs=np.arange(0,200,25)
            elif hvthr==50:
                if mmean == 'tot':
                    clevs=np.arange(0,2000,400)
                elif mmean =='mon':
                    clevs=np.arange(0,100,10)

        elif ptype=='aper_wet_cnt':
            cm = plt.cm.gnuplot2
            cbar_lab = '%'
            if hvthr==0:
                clevs=np.arange(0,90,10)
            elif hvthr==10:
                clevs=np.arange(0,45,5)
            elif hvthr==25:
                clevs=np.arange(0,10,1)
            elif hvthr==50:
                clevs=np.arange(0,4.5,0.5)

        elif ptype=='aper_wet_sum':
            cm = plt.cm.gnuplot2
            cbar_lab = '%'
            if hvthr==0:
                clevs=np.arange(0,100,10) # this should be 100% everywhere
            elif hvthr==10:
                clevs=np.arange(50,95,5)
            elif hvthr==25:
                clevs=np.arange(0,70,10)
            elif hvthr==50:
                clevs=np.arange(0,40,5)


        elif ptype=='tper_wet_cnt':
            cm = plt.cm.gnuplot2
            cbar_lab = '%'
            clevs = np.arange(0, 100, 10)

        elif ptype=='tper_wet_sum':
            cm = plt.cm.gnuplot2
            cbar_lab = '%'
            if hvthr==0:
                #clevs = [0, 10, 20, 30, 40, 50, 60] # should be the same as per_ttt
                clevs = np.arange(0, 60, 10)
            else:
                clevs=np.arange(0,60,10)

        elif ptype=='per_tttd_wet':
            cm = plt.cm.gnuplot2
            cbar_lab = '%'
            if hvthr==0:
                clevs = np.arange(0, 100, 10)
            elif hvthr==10:
                clevs = np.arange(0,50,5)
            elif hvthr==25:
                clevs = np.arange(0,30,3)
            elif hvthr==50:
                clevs = np.arange(0,15,1)


        # ttt cbars
        elif ptype=='tot_ttt':
            if mmean=='tot':
                clevs=[0,250,500,750,1000,1250,1500,1750,2000]
            elif mmean=='mon':
                clevs=[0,20,40,60,80,100,120,140]
            elif mmean=='day':
                clevs=[0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0]
            cm=plt.cm.YlGnBu
            cbar_lab='mm'
        elif ptype=='per_ttt':
            clevs=[0,10,20,30,40,50,60]
            cm = plt.cm.gnuplot2
            cbar_lab='%'
        elif ptype=='rain_per_ttt':
            clevs = [0.0,1.5,3.0,4.5,6.0,7.5,9.0,10.5,12.0]
            cm = plt.cm.YlGnBu
            cbar_lab='mm'
        elif ptype=='comp_anom_ttt' or ptype == 'comp_anom_ag':
            clevs = np.arange(-4, 4.5, 0.5)
            cm = plt.cm.seismic_r
            cbar_lab='mm'
        elif ptype=='comp_anom_cnt':
            clevs= np.arange(0,50,5)
            cm = plt.cm.magma_r
            cbar_lab='%'


        if test:
            cs = m.contourf(plon, plat, data4plot, cmap=cm, extend='both')
        else:
            cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')

        if labels:
            if ptype=='tot_ttt' or ptype=='comp_anom_ttt' or ptype == 'comp_anom_ag':
                tit=stats.mndict[mn]+': '+str(nttt_mon)+' TTT days '+str(int(round(float(nttt_mon)/float(nys))))+'/yr'
            elif ptype=='per_ttt' or ptype=='tper_wet_sum':
                tit=stats.mndict[mn]+': '+str(int(round((float(nttt_mon)/float(ndays_mon))*100.0)))+'% of days have TTTs'
            elif ptype=='all_wet_cnt' or ptype=='all_wet_sum':
                tit=stats.mndict[mn]+': mean = '+str(int(round(titstat)))+' per mon'
            else:
                tit = stats.mndict[mn]
        else:
            tit=stats.mndict[mn]
        plt.title(tit)

        if ptype == 'comp_anom_ag':
            if nttt_mon >= 1:
                hatch = m.contourf(plon, plat, mask_zeros, levels=[-1.0, 0.0, 1.0], hatches=["", '.'], alpha=0)


        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(data4plot)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.05,wspace=0.2,hspace=0.2)
    axcl=g.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs, cax=axcl)
    cbar.set_label(cbar_lab)

    if savefig:
        plt.savefig(figdir+'/Rainmap_'+ptype+'_'+file_suffix+'_'+under_of+'.png',dpi=150)

    return data4plot

def gridolrmap_season(s,eventkeys,olr,lat,lon,dtime,cl,season='coreseason',key='noaa-olr-0-0',\
                       ptype='ave_all',mmean='mon',under_of='dayof',figdir='test',file_suffix='test',\
                       savefig=False, test=True, labels=False,agthresh='perc_ag'):
    '''Produces subplots of ttt olr by month
    need to open the olr data with lon and lat and also the synop file
    see e.g. plot_ttt_olr_autothresh.py

    plot types
        ave_all - ave olr (for reference)
        ave_ttt - ave olr from TTTs
        comp_anom_ttt - ave_ttt as anom from monthly mean
        comp_anom_ag - as comp anom but with ag test
        comp_anom_cnt - % of TTT days with +ve or -ve anomalies
    under_of -> "dayof" is rain on day of TTTs, "under" is rain under TTTs
    '''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    yrs = np.unique(dtime[:,0])
    nys=len(yrs)

    # Get list of dates for these events
    edts = []
    ecnt=1
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[key]['mbt'][e.ixflags]
        #dts=e.trkdtimes
        for dt in range(len(dts)):
            if ecnt==1:
                edts.append(dts[dt])
            else:
                tmpdt=np.asarray(edts)
                # Check if it exists already
                ix = my.ixdtimes(tmpdt, [dts[dt][0]], \
                     [dts[dt][1]], [dts[dt][2]], [0])
                if len(ix)==0:
                    edts.append(dts[dt])
            ecnt+=1
    edts = np.asarray(edts)
    edts[:, 3] = 0
    print "Number of original TTT days found =  " + str(len(edts))

    # Get n lat and lon
    nlon=len(lon)
    nlat=len(lat)

    # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season == 'dryseason':mns = [4, 5, 6, 7, 8, 9]
    elif isinstance(season,list):
        mns=season
    print season
    print mns

    # Draw basemap
    m, f = blb.SAfrBasemap(lat,lon,drawstuff=True,prj='cyl',fno=1,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    cnt=1
    msklist=[]
    for mn in mns:
        print 'Plotting month '+str(mn)

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data
        #Ave OLR

        # First select the years you want - keeping this to allow option to select years later
        firstyear=yrs[0]
        lastyear=yrs[nys-1]
        olrdat=np.where((dtime[:,0]>=firstyear) & (dtime[:,0]<=lastyear))
        olrperiod=olr[olrdat,:,:]
        olrperiod=np.squeeze(olrperiod)
        newdates=dtime[olrdat]

        # Then select the month
        ix=np.where((newdates[:,1]==mn))
        olrmon=olrperiod[ix,:,:]
        datesmon=newdates[ix]
        ndays_mon=len(datesmon)
        olrmon=np.squeeze(olrmon)
        olrsum_all=np.nansum(olrmon,0)
        olrave_monthly=olrsum_all/nys
        olrave_daily=olrsum_all/ndays_mon



        if ptype=='ave_all':
            if mmean=='mon':
                data4plot=olrave_monthly
            elif mmean=='day':
                data4plot=olrave_daily

        else:

            #TTT olr
            if under_of=='dayof':
                ix2=np.where((edts[:,1]==mn))
                edatesmon=edts[ix2]

                indices = []
                for edt in range(len(edatesmon)):
                    ix = my.ixdtimes(datesmon, [edatesmon[edt][0]], \
                                 [edatesmon[edt][1]], [edatesmon[edt][2]], [0])
                    if len(ix)>=1:
                        indices.append(ix)
                if len(indices)>=2:
                    indices = np.squeeze(np.asarray(indices))
                else:
                    indices = indices
                nttt_mon=len(indices)
                print "Number of TTT days found in olr dataset for mon "+str(mn)+" =  " + str(nttt_mon)

                if nttt_mon==0:
                    olrave_ttt=np.zeros((nlat,nlon),dtype=np.float32)
                else:
                    olrsel=olrmon[indices,:,:]
                    if nttt_mon >= 2:
                        olrave_ttt=np.nanmean(olrsel,0)
                    else:
                        olrave_ttt=np.squeeze(olrsel)

                comp_anom = olrave_ttt - olrave_daily

                if ptype=='comp_anom_ag' or ptype=='comp_anom_cnt':
                    if nttt_mon >= 1:

                        anoms = np.zeros((nttt_mon, nlat, nlon), dtype=np.float32)
                        for day in range(nttt_mon):
                            this_anom = olrsel[day, :, :] - olrave_daily
                            anoms[day, :, :] = this_anom

                        if ptype=='comp_anom_ag':

                            anoms_signs = np.sign(anoms)
                            comp_signs = np.sign(comp_anom)

                            mask_zeros = np.zeros((nlat, nlon), dtype=np.float32)
                            for i in range(nlat):
                                for j in range(nlon):
                                    count = len(np.where(anoms_signs[:, i, j] == comp_signs[i, j])[0])
                                    perc = (float(count) / float(nttt_mon)) * 100
                                    if perc >= agthresh:
                                        mask_zeros[i, j] = 1
                                    else:
                                        mask_zeros[i, j] = 0

                        elif ptype=='comp_anom_cnt':

                            pos_pcent=np.zeros((nlat, nlon), dtype=np.float32)
                            neg_pcent=np.zeros((nlat, nlon), dtype=np.float32)
                            #zero_pcent=np.zeros((nlat, nlon), dtype=np.float32)


                            for i in range(nlat):
                                for j in range(nlon):
                                    count_p = len(np.where(anoms[:, i, j] > 0)[0])
                                    count_n = len(np.where(anoms[:, i, j] < 0)[0])
                                    #count_z = len(np.where(anoms[:, i, j] == 0)[0])


                                    perc_p = (float(count_p) / float(nttt_mon)) * 100
                                    perc_n = (float(count_n) / float(nttt_mon)) * 100
                                    #perc_z = (float(count_z) / float(nttt_mon)) * 100
                                    #print perc_z


                                    pos_pcent[i,j]=perc_p
                                    neg_pcent[i,j]=perc_n
                                    #zero_pcent[i,j]=perc_z

            if ptype=='ave_ttt':
                data4plot = olrave_ttt
            elif ptype == 'comp_anom_ttt':
                data4plot = comp_anom
            elif ptype == 'comp_anom_ag':
                #data4plot = comp_anom[::2, ::2]
                data4plot = comp_anom
            elif ptype == 'comp_anom_cnt':
                data4plot= neg_pcent
                #data4plot= pos_pcent

            # if ptype == 'comp_anom_ag':
            #     newlon = lon[::2]
            #     newlat = lat[::2]
            # else:
            #     newlon = lon
            #     newlat = lat

            newlon = lon
            newlat = lat


        #Plot
        plon,plat = np.meshgrid(newlon,newlat)

        if ptype=='ave_all':
            clevs = np.arange(200, 280, 10)
            cm = plt.cm.gray_r
        elif ptype=='ave_ttt':
            clevs = np.arange(200, 280, 10)
            cm=plt.cm.gray_r
        elif ptype=='comp_anom_ttt' or ptype == 'comp_anom_ag':
            clevs = np.arange(-12, 14, 2)
            cm = plt.cm.BrBG_r
        elif ptype == 'comp_anom_cnt':
            clevs= np.arange(30,75,5)
            cm = plt.cm.BrBG
            #clevs = [50, 55, 60, 65, 70]
            #cm = plt.cm.OrRd
            #cm = plt.cm.PuBu

        if test:
            cs = m.contourf(plon, plat, data4plot, cmap=cm, extend='both')
        else:
            cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')
        if labels:
            if ptype=='ave_ttt' or ptype=='comp_anom_ttt' or ptype == 'comp_anom_ag':
                tit=stats.mndict[mn]+': '+str(nttt_mon)+' TTT days '+str(int(round(float(nttt_mon)/float(nys))))+'/yr'
        else:
            tit=stats.mndict[mn]
        plt.title(tit)

        if ptype == 'comp_anom_ag':
            if nttt_mon >= 1:
                #mask_zeros=mask_zeros[::2,::2]
                hatch = m.contourf(plon, plat, mask_zeros, levels=[-1.0, 0.0, 1.0], hatches=["", '.'], alpha=0)

        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(data4plot)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.05,wspace=0.2,hspace=0.2)
    axcl=g.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs, cax=axcl)
    if ptype=='comp_anom_cnt':
        cbar.set_label('% neg')
    else:
        cbar.set_label('W/m^2')

    if savefig:
        plt.savefig(figdir+'/OLRmap_'+ptype+'_'+file_suffix+'_'+under_of+'.png',dpi=150)

    return data4plot

def gridvarmap_season(s,eventkeys,varstr,vardata,varlat,varlon,vardtime,olrunits,olrcl,season='coreseason',key='noaa-olr-0-0',\
                       ptype='comp_anom_ag',under_of='dayof',figdir='test',file_suffix='test',\
                       savefig=False, test=True, agthresh='perc_ag'):
    '''Produces composite plots assoc with TTTs
    need to open the var data with lon and lat and also the synop file
    see e.g. plot_ttt_omega_autothresh.py

    plot types
        comp_anom_ag - as comp anom but with ag test
        comp_anom_cnt - % of TTT days with +ve or -ve anomalies
    under_of -> "dayof" is rain on day of TTTs, "under" is rain under TTTs
    '''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    yrs = np.unique(vardtime[:,0])
    nys=len(yrs)

    # Get list of dates for these events
    edts = []
    ecnt=1
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[key]['mbt'][e.ixflags]
        #dts=e.trkdtimes
        for dt in range(len(dts)):
            if ecnt==1:
                edts.append(dts[dt])
            else:
                tmpdt=np.asarray(edts)
                # Check if it exists already
                ix = my.ixdtimes(tmpdt, [dts[dt][0]], \
                     [dts[dt][1]], [dts[dt][2]], [0])
                if len(ix)==0:
                    edts.append(dts[dt])
            ecnt+=1
    edts = np.asarray(edts)
    edts[:, 3] = 0
    print "Number of original TTT days found =  " + str(len(edts))

    # Get n lat and lon
    nlon=len(varlon)
    nlat=len(varlat)

    # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season == 'dryseason':mns = [4, 5, 6, 7, 8, 9]
    elif isinstance(season,list):
        mns=season
    print season
    print mns

    # Draw basemap
    m, f = blb.SAfrBasemap(varlat,varlon,drawstuff=True,prj='cyl',fno=1,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    cnt=1
    msklist=[]
    for mn in mns:
        print 'Plotting month '+str(mn)

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data

        # First select the years you want - keeping this to allow option to select years later
        firstyear=yrs[0]
        lastyear=yrs[nys-1]
        vardat=np.where((vardtime[:,0]>=firstyear) & (vardtime[:,0]<=lastyear))
        varperiod=vardata[vardat,:,:]
        varperiod=np.squeeze(varperiod)
        newdates=vardtime[vardat]

        # Then select the month
        ix=np.where((newdates[:,1]==mn))
        varmon=varperiod[ix,:,:]
        datesmon=newdates[ix]
        ndays_mon=len(datesmon)
        varmon=np.squeeze(varmon)
        varsum_all=np.nansum(varmon,0)
        varave_monthly=varsum_all/nys
        varave_daily=varsum_all/ndays_mon

        #Select TTT days
        if under_of=='dayof':
            ix2=np.where((edts[:,1]==mn))
            edatesmon=edts[ix2]

            indices = []
            for edt in range(len(edatesmon)):
                ix = my.ixdtimes(datesmon, [edatesmon[edt][0]], \
                             [edatesmon[edt][1]], [edatesmon[edt][2]], [0])
                if len(ix)>=1:
                    indices.append(ix)
            if len(indices)>=2:
                indices = np.squeeze(np.asarray(indices))
            else:
                indices = indices
            nttt_mon=len(indices)
            print "Number of TTT days found in dataset for mon "+str(mn)+" =  " + str(nttt_mon)

            if nttt_mon==0:
                varave_ttt=np.zeros((nlat,nlon),dtype=np.float32)
            else:
                varsel=varmon[indices,:,:]
                if nttt_mon >= 2:
                    varave_ttt=np.nanmean(varsel,0)
                else:
                    varave_ttt=np.squeeze(varsel)

            comp_anom = varave_ttt - varave_daily

            if ptype=='comp_anom_ag' or ptype=='comp_anom_cnt':
                if nttt_mon >= 1:

                    anoms = np.zeros((nttt_mon, nlat, nlon), dtype=np.float32)
                    for day in range(nttt_mon):
                        this_anom = varsel[day, :, :] - varave_daily
                        anoms[day, :, :] = this_anom

                    if ptype=='comp_anom_ag':

                        anoms_signs = np.sign(anoms)
                        comp_signs = np.sign(comp_anom)

                        mask_zeros = np.zeros((nlat, nlon), dtype=np.float32)
                        for i in range(nlat):
                            for j in range(nlon):
                                count = len(np.where(anoms_signs[:, i, j] == comp_signs[i, j])[0])
                                perc = (float(count) / float(nttt_mon)) * 100
                                if perc >= agthresh:
                                    mask_zeros[i, j] = 1
                                else:
                                    mask_zeros[i, j] = 0

                    elif ptype=='comp_anom_cnt':

                        pos_pcent=np.zeros((nlat, nlon), dtype=np.float32)
                        neg_pcent=np.zeros((nlat, nlon), dtype=np.float32)
                        #zero_pcent=np.zeros((nlat, nlon), dtype=np.float32)


                        for i in range(nlat):
                            for j in range(nlon):
                                count_p = len(np.where(anoms[:, i, j] > 0)[0])
                                count_n = len(np.where(anoms[:, i, j] < 0)[0])
                                #count_z = len(np.where(anoms[:, i, j] == 0)[0])


                                perc_p = (float(count_p) / float(nttt_mon)) * 100
                                perc_n = (float(count_n) / float(nttt_mon)) * 100
                                #perc_z = (float(count_z) / float(nttt_mon)) * 100
                                #print perc_z


                                pos_pcent[i,j]=perc_p
                                neg_pcent[i,j]=perc_n
                                #zero_pcent[i,j]=perc_z

            if ptype == 'comp_anom_ag':
                data4plot = comp_anom
            elif ptype == 'comp_anom_cnt':
                data4plot= neg_pcent

            newlon = varlon
            newlat = varlat


        #Plot
        plon,plat = np.meshgrid(newlon,newlat)

        if ptype=='comp_anom_ag':
            if varstr=='omega':
                clevs = np.arange(-0.06, 0.07, 0.01)
                cm = plt.cm.seismic
            elif varstr=='q':
                clevs = np.arange(-0.002, 0.00215, 0.00025)
                cm = plt.cm.seismic_r

        elif ptype == 'comp_anom_cnt':
            if varstr == 'omega':
                clevs= np.arange(30,75,5)
                cm = plt.cm.seismic_r
            elif varstr == 'q':
                clevs= np.arange(30,75,5)
                cm = plt.cm.seismic

        if test:
            cs = m.contourf(plon, plat, data4plot, cmap=cm, extend='both')
        else:
            cs = m.contourf(plon, plat, data4plot, clevs, cmap=cm, extend='both')
        tit=stats.mndict[mn]
        plt.title(tit)

        if ptype == 'comp_anom_ag':
            if nttt_mon >= 1:
                hatch = m.contourf(plon, plat, mask_zeros, levels=[-1.0, 0.0, 1.0], hatches=["", '.'], alpha=0)

        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(data4plot)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.05,wspace=0.2,hspace=0.2)
    axcl=g.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs, cax=axcl)
    if ptype=='comp_anom_cnt':
        cbar.set_label('% neg')

    if savefig:
        plt.savefig(figdir+'/Map_'+ptype+'_'+file_suffix+'_'+under_of+'.png',dpi=150)

    return data4plot

def gridvectmap_season(s,eventkeys,varstr,vardata_u,vardata_v,varlat,varlon,vardtime,olrunits,olrcl,season='coreseason',key='noaa-olr-0-0',\
                       ptype='comp_anom_ag',under_of='dayof',figdir='test',file_suffix='test',\
                       savefig=False, test=True, agthresh='perc_ag'):
    '''Produces composite plots assoc with TTTs
    need to open the var data with lon and lat and also the synop file
    see e.g. plot_ttt_qflux_autothresh.py

    plot types
        comp_anom_ttt - comp anom without ag test
        comp_anom_ag - as comp anom but with ag test
    under_of -> "dayof" is rain on day of TTTs, "under" is rain under TTTs
    '''
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    yrs = np.unique(vardtime[:,0])
    nys=len(yrs)

    # Get list of dates for these events
    edts = []
    ecnt=1
    for k in eventkeys:
        e = s.events[k]
        dts = s.blobs[key]['mbt'][e.ixflags]
        #dts=e.trkdtimes
        for dt in range(len(dts)):
            if ecnt==1:
                edts.append(dts[dt])
            else:
                tmpdt=np.asarray(edts)
                # Check if it exists already
                ix = my.ixdtimes(tmpdt, [dts[dt][0]], \
                     [dts[dt][1]], [dts[dt][2]], [0])
                if len(ix)==0:
                    edts.append(dts[dt])
            ecnt+=1
    edts = np.asarray(edts)
    edts[:, 3] = 0
    print "Number of original TTT days found =  " + str(len(edts))

    # Get n lat and lon
    nlon=len(varlon)
    nlat=len(varlat)

    # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
        elif season == 'dryseason':mns = [4, 5, 6, 7, 8, 9]
    elif isinstance(season,list):
        mns=season
    print season
    print mns

    # Draw basemap
    m, f = blb.SAfrBasemap(varlat,varlon,drawstuff=True,prj='cyl',fno=1,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,12])
    cnt=1
    msklist=[]
    for mn in mns:
        print 'Plotting month '+str(mn)

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data

        # First select the years you want - keeping this to allow option to select years later
        firstyear=yrs[0]
        lastyear=yrs[nys-1]
        vardat=np.where((vardtime[:,0]>=firstyear) & (vardtime[:,0]<=lastyear))

        varperiod_u=np.squeeze(vardata_u[vardat,:,:])
        varperiod_v=np.squeeze(vardata_v[vardat,:,:])

        newdates=vardtime[vardat]

        # Then select the month
        ix=np.where((newdates[:,1]==mn))

        varmon_u=np.squeeze(varperiod_u[ix,:,:])
        varmon_v=np.squeeze(varperiod_v[ix,:,:])

        datesmon=newdates[ix]
        ndays_mon=len(datesmon)

        varsum_all_u=np.nansum(varmon_u,0)
        varsum_all_v=np.nansum(varmon_v,0)


        varave_monthly_u=varsum_all_u/nys
        varave_monthly_v=varsum_all_v/nys

        varave_daily_u=varsum_all_u/ndays_mon
        varave_daily_v=varsum_all_v/ndays_mon

        #Select TTT days
        if under_of=='dayof':
            ix2=np.where((edts[:,1]==mn))
            edatesmon=edts[ix2]

            indices = []
            for edt in range(len(edatesmon)):
                ix = my.ixdtimes(datesmon, [edatesmon[edt][0]], \
                             [edatesmon[edt][1]], [edatesmon[edt][2]], [0])
                if len(ix)>=1:
                    indices.append(ix)
            if len(indices)>=2:
                indices = np.squeeze(np.asarray(indices))
            else:
                indices = indices
            nttt_mon=len(indices)
            print "Number of TTT days found in dataset for mon "+str(mn)+" =  " + str(nttt_mon)

            if nttt_mon==0:
                varave_ttt_u=np.zeros((nlat,nlon),dtype=np.float32)
                varave_ttt_v=np.zeros((nlat,nlon),dtype=np.float32)

            else:
                varsel_u=varmon_u[indices,:,:]
                varsel_v=varmon_v[indices,:,:]

                if nttt_mon >= 2:
                    varave_ttt_u=np.nanmean(varsel_u,0)
                    varave_ttt_v=np.nanmean(varsel_v,0)

                else:
                    varave_ttt_u=np.squeeze(varsel_u)
                    varave_ttt_v=np.squeeze(varsel_v)


            comp_anom_u = varave_ttt_u - varave_daily_u
            comp_anom_v = varave_ttt_v - varave_daily_v


            if ptype=='comp_anom_ag' or ptype=='comp_anom_cnt':
                if nttt_mon >= 1:

                    anoms_u = np.zeros((nttt_mon, nlat, nlon), dtype=np.float32)
                    anoms_v = np.zeros((nttt_mon, nlat, nlon), dtype=np.float32)

                    for day in range(nttt_mon):
                        this_anom_u = varsel_u[day, :, :] - varave_daily_u
                        this_anom_v = varsel_v[day, :, :] - varave_daily_v


                        anoms_u[day, :, :] = this_anom_u
                        anoms_v[day, :, :] = this_anom_v


                    if ptype=='comp_anom_ag':

                        anoms_signs_u = np.sign(anoms_u)
                        anoms_signs_v = np.sign(anoms_v)

                        comp_signs_u = np.sign(comp_anom_u)
                        comp_signs_v = np.sign(comp_anom_v)


                        mask_zeros = np.zeros((nlat, nlon), dtype=np.float32)
                        for i in range(nlat):
                            for j in range(nlon):
                                count_u = len(np.where(anoms_signs_u[:, i, j] == comp_signs_u[i, j])[0])
                                count_v = len(np.where(anoms_signs_v[:, i, j] == comp_signs_v[i, j])[0])

                                perc_u = (float(count_u) / float(nttt_mon)) * 100
                                perc_v = (float(count_v) / float(nttt_mon)) * 100

                                if perc_u >= agthresh and perc_v >= agthresh:
                                    mask_zeros[i, j] = 1
                                else:
                                    mask_zeros[i, j] = 0

                        # Set masked as 0
                        zeroed_comp_u = comp_anom_u * mask_zeros
                        zeroed_comp_v = comp_anom_v * mask_zeros



            if ptype =='comp_anom_ttt':
                data4plot_u = comp_anom_u
                data4plot_v = comp_anom_v


            elif ptype == 'comp_anom_ag':
                data4plot_u=zeroed_comp_u
                data4plot_v=zeroed_comp_v

            newlon = varlon
            newlat = varlat


        #Plot
        plon,plat = np.meshgrid(newlon,newlat)

        wind_sc = 1
        usc = 0.01
        lab = '0.01 kg/kg/ms'

        q = plt.quiver(newlon, newlat, data4plot_u, data4plot_v, scale=wind_sc)
        if cnt == 1:
            plt.quiverkey(q, X=0.9, Y=1.1, U=usc, label=lab, labelpos='W', fontproperties={'size': 'xx-small'})

        tit=stats.mndict[mn]
        plt.title(tit)

        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(data4plot_u)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.05,wspace=0.2,hspace=0.2)

    if savefig:
        plt.savefig(figdir+'/Map_'+ptype+'_'+file_suffix+'_'+under_of+'.png',dpi=150)

    return data4plot_u


def gridrainmap_bias_season(s,raingrid,rain,lat,lon,rdtime,eventkeys,yrs,figno=1,season='coreseason',key='um-olr-0-0',ptype='diff_tot',file_suffix='test',savefig=False):
    '''Produces subplots of ttt rainfall by month
    need to open the rain data with lon and lat and also the synop file
    see e.g. rainmap_test_months.py

    plot types
        diff_tot - UM total rain - TRMM total rain
        diff_ttt - UM ttt rain - TRMM ttt rain
        percent - (diff_ttt/diff_tot)*100 - % of bias from ttt bias

    '''
    mbkl=key.split('-')
    if not eventkeys:
        eventkeys=[]
        for ed in s.uniques:
            eventkeys.append(ed[0])

    if yrs == None:
        print "No years specified"
        yrs = np.unique(rdtime[:,0])
    else:
        print "years specified"
        print yrs
        yrs=yrs
    print yrs

    nlon=len(lon)
    nlat=len(lat)

    # Open TRMM data
    refdir='/ouce-home/students/bras2392/work/impala/data/refdata/'
    trmmfile=refdir+'trmm_3b42v7_n216.daily.precip.1998_2013.mon.nc'

    train, time, tlat, tlon, tdtime = \
        mync.opentrmm(trmmfile,'pcp','trmm',subs='SA_TR')
    tdtime[:,3]=0 # editing the hour time in date
    traingrid=(train,tdtime,(tlon,tlat))

    tnlon=len(tlon)
    tnlat=len(tlat)

    # Open satellite events
    satdir='/ouce-home/students/bras2392/PycharmProjects/ttcb_um_test/test/'
    satsyfile = satdir + "NOAA-OLR.synop"
    sats = sy.SynopticEvents((),[satsyfile],COL=False)
    satks = sats.events.keys();satks.sort()
    tkey='noaa-olr-0-0'
    tmbkl=tkey.split('-')

   # Get months
    if isinstance(season,str):
        if season=='coreseason':mns=[10,11,12,1,2,3]
        elif season=='fullseason':mns=[8,9,10,11,12,1,2,3,4,5,6,7]
    elif isinstance(season,list):
        mns=season

    # Draw basemap
    m, f = blb.SAfrBasemap(lat,lon,drawstuff=True,prj='cyl',fno=figno,rsltn='l')

    # Get multiplot
    if len(mns)==12:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    elif len(mns)==6:
        plt.close()
        g, axls = plt.subplots(figsize=[12,10])
    cnt=1
    msklist=[]
    for mn in mns:
        print mn

        if len(mns)==12:plt.subplot(4,3,cnt)
        elif len(mns)==6:plt.subplot(3,2,cnt)

        # Get the plot data (allmask)
        speckeys = stats.specificmon(s,eventkeys,yrs,mn,mbkl[0])
        tspeckeys = stats.specificmon(sats,satks,yrs,mn,tmbkl[0])

        #All rain
        firstyear=yrs[0]
        lastyear=yrs[len(yrs)-1]
        raindat=np.where((rdtime[:,0]>=firstyear) & (rdtime[:,0]<=lastyear))
        traindat=np.where((tdtime[:,0]>=firstyear) & (tdtime[:,0]<=lastyear))

        rainperiod=rain[raindat,:,:]
        trainperiod=train[traindat,:,:]

        rainperiod=np.squeeze(rainperiod)
        trainperiod=np.squeeze(trainperiod)

        newdates=rdtime[raindat]
        tnewdates=tdtime[traindat]

        ix=np.where((newdates[:,1]==mn))
        tix=np.where((tnewdates[:,1]==mn))

        rainmon=rainperiod[ix,:,:]
        trainmon=trainperiod[tix,:,:]

        rainmon=np.squeeze(rainmon)
        trainmon=np.squeeze(trainmon)

        rainsum_all=np.nansum(rainmon,0)
        trainsum_all=np.nansum(trainmon,0)

        allrainsum_diff=rainsum_all-trainsum_all

        if ptype=='diff_tot':
            allmask=allrainsum_diff
        else:
            #TTT rain
            rain4plot = stats.griddedrainmasks(s,speckeys,raingrid,refkey=key)
            train4plot = stats.griddedrainmasks(sats,tspeckeys,traingrid,refkey=tkey)

            nevent=rain4plot.shape[0]
            tnevent=train4plot.shape[0]

            data2sum=np.zeros((nevent,nlat,nlon),dtype=np.float32)
            tdata2sum=np.zeros((tnevent,tnlat,tnlon),dtype=np.float32)

            for e in xrange(nevent):
                data=rain4plot[e]
                ntstep=len(data[:,0,0])
                data2sum[e,:,:]=np.nansum(data,0)

            for e in xrange(tnevent):
                data=train4plot[e]
                ntstep=len(data[:,0,0])
                tdata2sum[e,:,:]=np.nansum(data,0)

            rainsum_ttt=np.nansum(data2sum,0)
            trainsum_ttt=np.nansum(tdata2sum,0)

            tttrainsum_diff=rainsum_ttt-trainsum_ttt

            if ptype=='diff_ttt':
                allmask=tttrainsum_diff
            elif ptype=='percent':
                biasperc_ttt=(tttrainsum_diff/allrainsum_diff)*100
                allmask=biasperc_ttt
            elif ptype=='rainperttt':
                dclim=(1,9,1)
                spatiofreq=stats.spatiofreq2(m,s,lat,lon,yrs,speckeys,clim=dclim,key=key,month=mn,flagonly=True)

        #Plot
        plon,plat = np.meshgrid(lon,lat)
        if ptype=='diff_tot':
            clevs=[-1500,-1000,-500,0,500,1000,1500]
            #cm=plt.cm.RdBu
            cm=plt.cm.cubehelix_r
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        elif ptype=='diff_ttt':
            clevs=[-450,-300,-150,0,150,300,450]
            cm=plt.cm.RdBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='both')
        elif ptype=='percent':
            clevs=[0,10,20,30,40,50,60,70,80,90,100]
            #Mask negative values
            allmask=np.ma.masked_where(allmask < 0, allmask)
            cm=plt.cm.RdPu
            cm.set_bad(color='white')
            cs = m.contourf(plon,plat,allmask,clevs,cmap=cm,extend='max')
        else:
            cm=plt.cm.RdYlBu
            #cm=plt.cm.cubehelix
            cs = m.contourf(plon,plat,allmask,cmap=cm)
            plt.colorbar()
        plt.title(monthname[mn-1])

        # redraw - only label latitudes if plot is on left
        if len(mns)==12:
            if cnt == 1 or cnt == 4 or cnt == 7 or cnt == 10:
                syp.redrawmap(m,lns=True,resol='verylow')
            else:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
        elif len(mns)==6:
            if cnt%2==0:
                syp.redrawmap(m,lns=True,resol='verylow',parallel=False)
            else:
                syp.redrawmap(m,lns=True,resol='verylow')
        cnt+=1
        msklist.append(allmask)
    plt.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.02,wspace=0.02,hspace=0.2)
    axcl=g.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs,cax=axcl)
    if ptype=='percent':cbar.set_label('%')
    else:cbar.set_label('mm')

    if savefig:
        plt.savefig('Grid_bias_TTTrain_months_'+file_suffix+'.png',dpi=150)

    return msklist
