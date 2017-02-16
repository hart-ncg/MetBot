# mynetcdf.py
'''Some useful netcdf file reading/writing routines'''

import numpy as np
try:
    import Scientific.IO.NetCDF as kh
    from mpl_toolkits.basemap import date2index, date2num, num2date
    print "Importing Scientific.IO.NetCDF library for netcdf."
    print "mpl_toolkits.basemap used for nc date support."
except ImportError:
    print "Importing netCDF4 library for netcdf & date support."
    import netCDF4 as kh
    kh.NetCDFFile = kh.Dataset
    date2index,date2num,num2date = kh.date2index,kh.date2num,kh.num2date
import time as tm
import datetime
import MetBot.dset_dict as dsetdict

# Add to this dictionary as need be by looking at ncdump -h ????.nc
dimdict={"ncep2": ['time','lat','lon','level','lev'],
"ncep": ['time','lat','lon'],
"20cr": ['time','lat','lon'],
"interp_olr": ['time','lat','lon'],
"noaa": ['time','lat','lon'],
"tamsat_ccd": ['time','lat','lon'],
"had": ['time','latitude','longitude','level'],
"um": ['t','latitude','longitude','toa'],
"umpr": ['t','latitude','longitude','surface'],
"cmip5": ['time','lat','lon'],
#"era": ['date','latitude','longitude','level'],
"era": ['time','latitude','longitude'],
"hadam3p": ['t','latitude','longitude','p','theta','surface','toa'],
"cfsr": ['time','latitude','longitude','level'],
"nddiagnc": ['t','y','x','p'],
"trmm": ['time','latitude','longitude']}
writedict={}

def isubs(sub,lat,lon,*args):
    ''' (ilat1, ilat2), (ilon1, ilon2), (ilev)  = isub(sub,lat,lon [,lev, levselect])

    Simple function to subset data

    INPUT: lat (1-dim array)
            lon (1-dim array)
            lev, levselect (optional)

            sub ((latmin,latmax),(lonmin,lonmax)) tuple
            or
            sub can be "Fauch08", "SH", "imkf", "playing
            or
            sub can be False: this behaviour will return (0,-1)
                              as indices for lat,lon but means
                              if isubs call included [lev,levsel] 
                              an level select indice will be return too"

    RETURNS: subdata, sublat, sublon'''
    domains={}
    domains['Fauch08'] = ((-40.0,-15.0),(7.5,70.0),)
    domains['Fauch08_mod'] = ((-40.0,-10.0),(7.5,70.0))
    domains['SH'] = ((-90.0,0.0),(0,360))
    domains['SH60'] = ((-60.0,0.0),(0,360))
    domains['imkf'] = ((-46.25,-10.0),(0.0,91.875))
    domains['SA'] = ((-60.0,0.0),(0.0,100.0))
    domains['SA_TR'] = ((-40.0,-15.0),(7.5,100.0))
    domains['NA'] = ((0.0,60.0),(280.0,359.0))
    domains['NP'] = ((0.0,60.0),(180.0,260.0))
    domains['SAROI'] = ((-35,-23),(7.5, 80))
    domains['WPR'] = ((-40.0,-15.0),(7.5, 40.0))
    domains['EPR'] = ((-40.0,-15.0),(40.0, 100.0))
    if isinstance(sub,str):
        domain=domains[sub]; getisubs=True
    elif isinstance(sub,tuple):
        domain=sub; getisubs=True
    elif isinstance(sub,bool):
        ilat12,ilon12=-99,-99 # due to indexing style in python -2 is
        # placed here so that when 1 is added (happens within netcdf.opennc)
        # indices are 0,-1 i.e. everything
        getisubs=False
    else:
        print "ERROR: sub needs to be a tuple or a string or a bool"
    if getisubs:
        lt1,lt2 = domain[0]; ln1,ln2 = domain[1]
        ilats = ((lt1<=lat) & (lat<=lt2)).nonzero()[0]
        ilons = ((ln1<=lon) & (lon<=ln2)).nonzero()[0]
        ilat12=(ilats[0],ilats[-1]);ilon12=(ilons[0],ilons[-1])
        ilev12=()
    if len(args)>0:
        lev, levselect = args
        if isinstance(levselect,tuple):
            if lev[-1]>lev[0]: lv1,lv2=levselect    # for level dim increasing 
            elif lev[0]>lev[-1]: lv2,lv1=levselect  # for level dim decreasing
            ilev = ((lv1<=lev) & (lev<=lv2)).nonzero()[0]
            ilev12=(ilev[0],ilev[-1])
        else:
            ilev=np.where(np.absolute((lev-levselect)).min() ==\
                           lev-levselect)[0][0]
            ilev12=(ilev,ilev)
        #print "Closest match for level:",lev[ilev]

    return ilat12, ilon12, ilev12

def ixdtimes(t,yrs,mns,ds,hs,mode='range'):
    '''ixt_1, ixt_2 = ixdtimes(t,yrs,mns,ds,hs)'''
    if mode=='match':
        teststr='';ytest=[];mtest=[];dtest=[];htest=[]
        if yrs: 
            ytest=["(t[:,0]=="+str(y)+")" for y in yrs]
            teststr=" | ".join(ytest) + " & "
        if mns: 
            mtest=["(t[:,1]=="+str(m)+")" for m in mns]
            teststr=teststr + " | ".join(mtest) + " & "
        if ds:
            dtest=["(t[:,2]=="+str(d)+")" for d in ds]
            teststr=teststr + " | ".join(dtest) + " & "
        if hs:
            htest=["(t[:,3]=="+str(h)+")" for h in hs]
            teststr=teststr + " | ".join(htest)
        if teststr[-2:-1] == "&":
            teststr=teststr[:-2]
        print 'ix = np.where('+teststr+')[0]'
        exec('ix = np.where('+teststr+')[0]')

    if mode=='range':
        teststr='';ytest=[];mtest=[];dtest=[];htest=[]
        tmin=yrs[0]*10e6 + mns[0]*10e4 + ds[0]*10e2 + hs[0];print tmin
        tmax=yrs[1]*10e6 + mns[1]*10e4 + ds[1]*10e2 + hs[1];print tmax

        tnum=t[:,0]*10e6 + t[:,1]*10e4 + t[:,2]*10e2 + t[:,3]
        ix = np.where((tnum >= tmin) & (tnum <= tmax))[0]

    return ix[0], ix[-1]+1

def fix360d(dtime):
    '''Quick fix for a problem seemed to be having with the 360-day
    hadam3p files'''
    for i in xrange(len(dtime)):
        if dtime[i].month==13:
            dtime[i].month=01
            #print "month changed to 01 since next on is, ",dtime[i+1]
            #try:
                    #print dtime[i]
            #except:
                    #print "Not fixed"
    return dtime

def dtime2arr(dtime):
    '''humandates = dtime2arr(dtime)
    
    Nifty tool converting arrays of datetime objects into arrays of 
    entries [YYYY,MM,DD,HH] for easy searching, index, etc.'''
    ymdh=np.zeros((len(dtime),4),dtype=np.int)
    for t in xrange(len(dtime)):
        ymdh[t,:] = dtime[t].timetuple()[:4]

    return ymdh

def addscalefill(ncf,data,varstr):
    '''This will do nothing if the netcdf module is netCDF4 as that
       much of this is built in there.'''
    if kh.__name__ == 'netCDF4':
        if isinstance(data,np.ma.MaskedArray):
            data = np.where(data.mask,np.nan,data.data)
        return data
    attrbs=[]
    exec('attrbs.append(getattr(ncf.variables[\''+varstr+'\']\
         ,\'missing_value\',False))')
    exec('attrbs.append(getattr(ncf.variables[\''+varstr+'\']\
         ,\'_FillValue\',False))')
    exec('attrbs.append(getattr(ncf.variables[\''+varstr+'\']\
         ,\'scale_factor\',False))')
    exec('attrbs.append(getattr(ncf.variables[\''+varstr+'\']\
         ,\'add_offset\',False))')
    if not attrbs[2]: attrbs[2]=1; #print "No scale factor"
    if not attrbs[3]: attrbs[3]=0; #print "No add_offset"
    # Replace missing values (attrbs[0]) or _FillValue (attrbs[1])
    if attrbs[0]:
        try:
            data = np.where(data == attrbs[0][0],np.NaN,data)
        except MemoryError:
            print "File too big, looping through time..."
            for it in xrange(data.shape[0]):
                data[it,:]=np.where(data[it,:]==attrbs[0][0],np.NaN,data[it,:])
        except IndexError:
            print "Nothing to NaN"
    elif attrbs[1]:
        try:
            data = np.where(data == attrbs[1][0],np.NaN,data)
        except MemoryError:
            print "File too big, looping through time..."
            for it in xrange(data.shape[0]):
                data[it,:]=np.where(data[it,:]==attrbs[1][0],np.NaN,data[it,:])
    # Apply scale factors
    try:
        data = data*attrbs[2] + attrbs[3]
    except MemoryError:
        for it in xrange(data.shape[0]):
            data[it,:] = data[it,:]*attrbs[2] + attrbs[3]

    return data


def opennc(ncfile,varstr,dset,sub=False,levselect=False,subtime=False):
    '''var, time, lat, lon, [[,lev], [,time_bnds]]= opennc(ncfile,varstr,dset,
                                        sub=False,levselect=False,subtime=False)

    Main function to open some observational and reanalysis data sets with ease.
    Flexibility comes from dimdict dictionary which holds dimension names,
    Option to open only subset lat, lon & lev which is necessary for 
    some large files.

    NOTE:- COARDS & CF compliance assumed variable.shape = time, [lev,] lat, lon
         - flexibility and exceptions for subsetting not fully implemented,
           will work but you may find buginess
         - hadam3p is .nc output from subset.tcl (uses xconv), where input was
           UM .pp files

    USAGE: varstr - string
           dset   - string valid: ncep2, interp_olr, had, era, hadam3p, cfsr
           sub    - tuple - ((latmin,latmax),(lonmin,lonmax))
                   or string - see options availble in mynetcdf.isubs
           levselect - will return level slice closest to given value
                       NOTE: to use this option without subsetting lat,lon grid
                             set sub=True
           subtime - tuple (yyyy, mm ,dd ,hh, ss) each entry can be a list of 
                    years to select

    RETURNS: var, lat, lon, lev'''

    dimlist = dimdict[dset][:]
    timestr = dimlist[0]
    ncf = kh.NetCDFFile(ncfile,'r')
    #vkeys = ncf.variables.keys() # could use something with this?
    for i in dimlist[:]:
        try:
            exec(i + ' = np.float32(ncf.variables[\''+ i + '\'][:])')
        except:
            print 'Variable \"'+i+ '\" does not exist in '+ncfile
            dimlist.remove(i);continue

    # HUMAN TIME CONVERSION AND TIME SUBSET IF REQUIRED \
    # (only really used for big files like cfsr)

    if dset=='hadam3p' or ncfile.split('/').count('flavours'):
        exec('dtime=num2date(('+timestr+'-1)*24,\
                 units="hours since 1959-12-01 00:00:00",calendar="360_day")') 
        # the above t-1 thing is a hack, but it works apparently for me
        dtime=fix360d(dtime)
    elif dset=='cfsr' and \
    ncf.variables['time'].units== "seconds since 1970-01-01 00:00:00.0 0:00":
        print "Performing hrtime conversion to hours since 1800-01-01"
        time=time/3600+1490184.0
        exec('dtime=num2date('+timestr+',\
             units="hours since 1800-01-01 00:00:00.0",calendar=\'gregorian\')')
    elif varstr=='precipitation' and not 'time' in locals():
        print "Giving time to timeless TRMM 3hr files..."
        day, hr =ncfile.split('/')[-1].split('.')[1].split('_')
        yyyy,mm,dd = day[:4],day[4:6],day[6:8]
        dtm = datetime.datetime(int(yyyy),int(mm),int(dd),int(hr))
        hrtm = date2num(dtm,units="hours since 1800-01-01 00:00:00.0",\
                        calendar='gregorian')
        dtime, time = np.asarray([dtm]), np.asarray([hrtm])
        dimlist.insert(0,'time')
    elif dset=='tamsat_ccd':
        fdate = ncfile.split('/')[-1]
        yr,mn,dy,thresh = fdate.split('_')
        yr=yr[3:]
        dtm = datetime.datetime(int(yr),int(mn),int(dy),0)
        hrtm = date2num(dtm,units="hours since 1800-01-01 00:00:00.0",\
                        calendar='gregorian')
        dtime, time = np.asarray([dtm]), np.asarray([hrtm])
    else:
        exec('dtime=num2date('+timestr+',units=ncf.variables[timestr].units,\
             calendar=\'gregorian\')')
    dtarr=dtime2arr(dtime)

    if not subtime:
        if sub and levselect:
            #print "lat, lon, lev subset"
            exec('ilats, ilons, ilev = isubs(sub,'+dimlist[1]+','+dimlist[2]+\
                 ','+dimlist[3]+',levselect)')
            if ilats==-99:
                ilt1, ilt2, iln1,iln2, ilv1, ilv2 = '0', '','0', '',\
                 str(ilev[0]),str(ilev[1]+1)
            else:
                ilt1,ilt2,iln1,iln2,ilv1,ilv2 = str(ilats[0]), str(ilats[1]+1),\
                  str(ilons[0]), str(ilons[1]+1),str(ilev[0]),str(ilev[1]+1) 
            exec('data = ncf.variables[\''+ varstr + '\'][:,'+ilv1+':'+ilv2+','\
                                      +ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
            exec(dimlist[3]+'='+dimlist[3]+'['+ilv1+':'+ilv2+']')
        elif sub:
            print "lat, lon subset"
            exec('ilats,ilons,ilev = isubs(sub,'+dimlist[1]+','+dimlist[2]+')')
            ilt1, ilt2, iln1,iln2 = str(ilats[0]), str(ilats[1]+1),\
                                    str(ilons[0]), str(ilons[1]+1)
            if len(dimlist)==3:
                if varstr=='precipitation':
                    exec('data = ncf.variables[\''+ varstr + '\']\
                                 ['+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
                else:
                    exec('data = ncf.variables[\''+ varstr + '\']\
                                 [:,'+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            elif len(dimlist)==4:
                exec('data = ncf.variables[\''+ varstr + '\']\
                             [:,:,'+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
        else:
            #print "No subsetting"
            exec('data = ncf.variables[\''+ varstr + '\'][:]')
    # DO SOME TIME SUBSETTING
    elif subtime:
        print "Opening subset of time..."
        yrs, mns, ds, hs = subtime
        it1, it2 = ixdtimes(dtarr,yrs,mns,ds,hs)
        dtarr=dtarr[it1:it2,:]
        it1, it2 = str(it1), str(it2)
        exec(timestr+'='+timestr+'['+it1+':'+it2+']')
        if sub and levselect:
            print "lat, lon, lev subset"
            exec('ilats, ilons, ilev = isubs(sub,'+dimlist[1]+','\
                  +dimlist[2]+','+dimlist[3]+',levselect)')
            if ilats==-99:
                ilt1, ilt2, iln1,iln2, ilv1, ilv2 = '0', '','0', '',\
                 str(ilev[0]),str(ilev[1]+1)
            else:
                ilt1,ilt2,iln1,iln2,ilv1,ilv2 = str(ilats[0]), str(ilats[1]+1),\
                  str(ilons[0]), str(ilons[1]+1),str(ilev[0]),str(ilev[1]+1) 
            exec('data = ncf.variables[\''+ varstr + '\'][:,'+ilv1+':'+ilv2+','\
                                      +ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
            exec(dimlist[3]+'='+dimlist[3]+'['+ilv1+':'+ilv2+']')
        elif sub:
            print "lat, lon subset"
            exec('ilats,ilons,ilev = isubs(sub,'+dimlist[1]+','+dimlist[2]+')')
            ilt1, ilt2, iln1,iln2 = str(ilats[0]), str(ilats[1]+1),\
                                        str(ilons[0]), str(ilons[1]+1)
            if len(dimlist)==3:
                exec('data = ncf.variables[\''+ varstr + '\']\
                         ['+it1+':'+it2+','+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            elif len(dimlist)==4:
                exec('data = ncf.variables[\''+ varstr + '\']\
                       ['+it1+':'+it2+',:,'+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
        else:
            #print "No subsetting"
            exec('data = ncf.variables[\''+ varstr + '\']['+it1+':'+it2+',:]')

    exec(varstr+' = addscalefill(ncf,data,varstr)')

    ncf.close()
    # Silly dataset specific tweaks
    if dset=='cfsr':
        latitude=latitude[::-1]
        exec(varstr+'='+varstr+'[:,::-1,:]')


    exec('out=(np.float32('+varstr+'),'+', '.join(dimlist)+',dtarr)')

    return out

def opennc2(ncfile,varstr,mname,dset,sub=False,levselect=False,subtime=False):
    '''var, time, lat, lon, [[,lev], [,time_bnds]]= opennc2(ncfile,varstr,mname,dset,
                                        sub=False,levselect=False,subtime=False)

    Main function to open some observational and reanalysis data sets with ease.
    Flexibility comes from dimdict dictionary which holds dimension names,
    Option to open only subset lat, lon & lev which is necessary for
    some large files.

    opennc2 is an edited version for greater flexibility for multiple models using dset_dict.py

    NOTE:- COARDS & CF compliance assumed variable.shape = time, [lev,] lat, lon
         - flexibility and exceptions for subsetting not fully implemented,
           will work but you may find buginess
         - hadam3p is .nc output from subset.tcl (uses xconv), where input was
           UM .pp files

    USAGE: varstr - string
           dset   - string valid: um, umpr, cmip5, noaa
           sub    - tuple - ((latmin,latmax),(lonmin,lonmax))
                   or string - see options availble in mynetcdf.isubs
           levselect - will return level slice closest to given value
                       NOTE: to use this option without subsetting lat,lon grid
                             set sub=True
           subtime - tuple (yyyy, mm ,dd ,hh, ss) each entry can be a list of
                    years to select

    RETURNS: var, lat, lon, lev'''

    dimlist = dimdict[dset][:]
    timestr = dimlist[0]
    ncf = kh.NetCDFFile(ncfile,'r')
    #vkeys = ncf.variables.keys() # could use something with this?
    for i in dimlist[:]:
        try:
            exec(i + ' = np.float32(ncf.variables[\''+ i + '\'][:])')
        except:
            print 'Variable \"'+i+ '\" does not exist in '+ncfile
            dimlist.remove(i);continue

    # HUMAN TIME CONVERSION AND TIME SUBSET IF REQUIRED \
    # (only really used for big files like cfsr)

    if dset=='um':
        moddct = dsetdict.dset_deets[dset][mname]
        units = moddct['timeunit']
        cal = moddct['calendar']
        if mname=='anqjn' or mname=='antib':
            exec ('dtime=num2date((' + timestr + '-1)*24,units="' + units + '",calendar="' + cal + '")')
        elif mname=='u-ab674' or mname=='u-ab680':
            exec('dtime=num2date((' + timestr + '-1),units="' + units + '",calendar="' + cal + '")')
        else:
            print 'new mname - check time string'
        # the above t-1 thing is a hack, but it works apparently for me
        dtime = fix360d(dtime)
    elif dset=='umpr':
        moddct = dsetdict.dset_deets[dset][mname]
        units = moddct['timeunit']
        cal = moddct['calendar']
        exec ('dtime=num2date((' + timestr + '-1)*24,units="' + units + '",calendar="' + cal + '")')
        # the above t-1 thing is a hack, but it works apparently for me
        dtime = fix360d(dtime)
    elif dset == 'cmip5':
        moddct = dsetdict.dset_deets[dset][mname]
        cal = moddct['calendar']
        units = moddct['timeunit']
        newunits = units.replace('days since', 'hours since')
        exec ('dtime=num2date((' + timestr + ')*24,units="' + newunits + '",calendar="' + cal + '")')
        if cal == '360_day':
            dtime = fix360d(dtime)
    else:
        exec('dtime=num2date('+timestr+',units=ncf.variables[timestr].units,\
             calendar=\'gregorian\')')
    dtarr=dtime2arr(dtime)

    if not subtime:
        if sub and levselect:
            #print "lat, lon, lev subset"
            exec('ilats, ilons, ilev = isubs(sub,'+dimlist[1]+','+dimlist[2]+\
                 ','+dimlist[3]+',levselect)')
            if ilats==-99:
                ilt1, ilt2, iln1,iln2, ilv1, ilv2 = '0', '','0', '',\
                 str(ilev[0]),str(ilev[1]+1)
            else:
                ilt1,ilt2,iln1,iln2,ilv1,ilv2 = str(ilats[0]), str(ilats[1]+1),\
                  str(ilons[0]), str(ilons[1]+1),str(ilev[0]),str(ilev[1]+1)
            exec('data = ncf.variables[\''+ varstr + '\'][:,'+ilv1+':'+ilv2+','\
                                      +ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
            exec(dimlist[3]+'='+dimlist[3]+'['+ilv1+':'+ilv2+']')
        elif sub:
            print "lat, lon subset"
            exec('ilats,ilons,ilev = isubs(sub,'+dimlist[1]+','+dimlist[2]+')')
            ilt1, ilt2, iln1,iln2 = str(ilats[0]), str(ilats[1]+1),\
                                    str(ilons[0]), str(ilons[1]+1)
            if len(dimlist)==3:
                if varstr=='precipitation':
                    exec('data = ncf.variables[\''+ varstr + '\']\
                                 ['+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
                else:
                    exec('data = ncf.variables[\''+ varstr + '\']\
                                 [:,'+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            elif len(dimlist)==4:
                exec('data = ncf.variables[\''+ varstr + '\']\
                             [:,:,'+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
        else:
            #print "No subsetting"
            exec('data = ncf.variables[\''+ varstr + '\'][:]')
    # DO SOME TIME SUBSETTING
    elif subtime:
        print "Opening subset of time..."
        yrs, mns, ds, hs = subtime
        it1, it2 = ixdtimes(dtarr,yrs,mns,ds,hs)
        dtarr=dtarr[it1:it2,:]
        it1, it2 = str(it1), str(it2)
        exec(timestr+'='+timestr+'['+it1+':'+it2+']')
        if sub and levselect:
            print "lat, lon, lev subset"
            exec('ilats, ilons, ilev = isubs(sub,'+dimlist[1]+','\
                  +dimlist[2]+','+dimlist[3]+',levselect)')
            if ilats==-99:
                ilt1, ilt2, iln1,iln2, ilv1, ilv2 = '0', '','0', '',\
                 str(ilev[0]),str(ilev[1]+1)
            else:
                ilt1,ilt2,iln1,iln2,ilv1,ilv2 = str(ilats[0]), str(ilats[1]+1),\
                  str(ilons[0]), str(ilons[1]+1),str(ilev[0]),str(ilev[1]+1)
            exec('data = ncf.variables[\''+ varstr + '\'][:,'+ilv1+':'+ilv2+','\
                                      +ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
            exec(dimlist[3]+'='+dimlist[3]+'['+ilv1+':'+ilv2+']')
        elif sub:
            print "lat, lon subset"
            exec('ilats,ilons,ilev = isubs(sub,'+dimlist[1]+','+dimlist[2]+')')
            ilt1, ilt2, iln1,iln2 = str(ilats[0]), str(ilats[1]+1),\
                                        str(ilons[0]), str(ilons[1]+1)
            if len(dimlist)==3:
                exec('data = ncf.variables[\''+ varstr + '\']\
                         ['+it1+':'+it2+','+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            elif len(dimlist)==4:
                exec('data = ncf.variables[\''+ varstr + '\']\
                       ['+it1+':'+it2+',:,'+ilt1+':'+ilt2+','+iln1+':'+iln2+']')
            exec(dimlist[1]+'='+dimlist[1]+'['+ilt1+':'+ilt2+']')
            exec(dimlist[2]+'='+dimlist[2]+'['+iln1+':'+iln2+']')
        else:
            #print "No subsetting"
            exec('data = ncf.variables[\''+ varstr + '\']['+it1+':'+it2+',:]')

    exec(varstr+' = addscalefill(ncf,data,varstr)')

    ncf.close()
    # Silly dataset specific tweaks
    if dset == 'cmip5':
        lat = lat[::-1]
        exec (varstr + '=' + varstr + '[:,::-1,:]')
    if dset == 'um':
        latitude=latitude[::-1]
        exec (varstr + '=' + varstr + '[:,:,::-1,:]')

    exec('out=(np.float32('+varstr+'),'+', '.join(dimlist)+',dtarr)')

    return out


def openncep2(ncfile,varstr,dset='ncep2',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def opennc_generic(ncfile,varstr,dset='ncep2',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def opentrmm(ncfile,varstr,dset='trmm',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def openhad(ncfile,varstr,dset='had',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def openera(ncfile,varstr,dset='era',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def openhadam3p(ncfile,varstr,dset='hadam3p',subs=False,levsel=False):
    vardct={'uwnd': 'u','vwnd':'v','hgt':'ht','olr':'olr','rhum':'rh','air':'t'}
    out = opennc(ncfile,vardct[varstr],dset,sub=subs,levselect=levsel)
    return out

def opencfsr(ncfile,varstr,dset='cfsr',subs=False,levsel=False,subt=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel,subtime=subt)
    return out

def openolr(ncfile,varstr,dset='interp_olr',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def opennddiagnc(ncfile,varstr,dset='nddiagnc',subs=False,levsel=False):
    out = opennc(ncfile,varstr,dset,sub=subs,levselect=levsel)
    return out

def opentamsatCCD(ncfile,varstr='ccd',dset='tamsat_ccd'):
    out = opennc(ncfile,varstr,dset)
    return out

def openolr_multi(ncfile,varstr,name,dataset='noaa',subs=False,levsel=False):
    out = opennc2(ncfile,varstr,name,dataset,sub=subs,levselect=levsel)
    return out


def openpscdf(ncfile,varstr,dset='pcdf',subs=False,levsel=False):
    psdimdict={"pcdf": ['time','dimy_T','dimx_T','dimz_T'],
    "scdf": ['time','dimy_PS','dimx_PS','dimz_TH']}
    dimlist = psdimdict[dset]
    timestr = dimlist[0]
    ncf = kh.NetCDFFile(ncfile,'r')
    #vkeys = ncf.variables.keys() # could use something with this?
    if kh.__name__=='Scientific.IO.NetCDF':
        for i in dimlist[:]:
            exec(i[2:4]+'=ncf.dimensions["'+i+'"]')
    else:
        for i in dimlist[:]:
            exec(i[2:4]+'=ncf.dimensions["'+i+'"].__len__()')
    varobj=ncf.variables[varstr]
    lon, lat, lev = np.linspace(varobj.xmin,varobj.xmax,mx), \
                    np.linspace(varobj.ymin,varobj.ymax,my),\
                    np.linspace(varobj.zmin,varobj.zmax,mz)
    lon=lon-360
    if levsel: 
        if isinstance(levsel,tuple):
            lv1,lv2=levsel
            ilev = ((lv2<=lev) & (lev<=lv1)).nonzero()[0]
            ilev12=(ilev[0],ilev[-1])
        else:
            ilev=np.where(np.absolute((lev-levsel)).min() ==\
                           lev-levsel)[0][0]
            ilev12=(ilev,ilev)
        ilv1,ilv2=ilev12
        vrb=varobj[0,ilv1:ilv2+1,:,:].squeeze()
        lev=lev[ilv1:ilv2+1]
    else: vrb=varobj[:].squeeze()
    ### Calculate true time
    day, hr =ncfile.split('/')[-1].split('_')
    yyyy,mm,dd = day[1:5],day[5:7],day[7:9]
    dtm = datetime.datetime(int(yyyy),int(mm),int(dd),int(hr))
    hrtm = date2num(dtm,units="hours since 1800-01-01 00:00:00.0",\
                    calendar='gregorian')
    dtime, time = np.asarray([dtm]), np.asarray([hrtm])

    out=(vrb, time, lat, lon, lev, dtime)
    ncf.close()
    return out
### WRITE OUT FILES
def writenc(ncfile,datastr,timeunits,dset,data,*args,**kwargs):
    '''writencep2(ncfile,datastr,timeunits,data,[time,lat,lon,level])

       This function writes out a NetCDF file with number of attributes
       included'''
    dimsnames=dimdict[dset]
    cntrl = len(args)
    outf=kh.NetCDFFile(ncfile,'w')
    data=np.float32(data)
    #outf.title = "EP type sensitivity experiment SST forcing field"
    outf.history = "Created with python module netCDF4 on " + tm.ctime()
    for ky in kwargs.keys():
        exec('outf.%s = "%s"' %(ky,kwargs[ky]))
    outf.createDimension('time',None)
    timef = outf.createVariable('time','f',('time',))
    #timef.units = "days since 1800-1-1 00:00:00"
    if timeunits:
        timef.units = timeunits
    else:
        timef.units = "hours since 1800-1-1 00:00:0.0"
    timef.long_name = "Time"
    timef.actual_range = [args[0][0],args[0][-1]]
    #timef.delta_t = "0000-01-00 00:00:00"
    #timef.prev_avg_period = "0000-00-07 00:00:00"
    #timef.standard_name = "time"
    timef.axis = "t"
    timef[:]=np.float32(args[0])

    if cntrl >= 2:
        outf.createDimension('lat',args[1].shape[0])
        latf = outf.createVariable('lat','f',('lat',))
        latf.units = "degrees_north"
        latf.long_name = "Latitude"
        latf.actual_range = [args[1][0], args[1][-1]]
        #latf.standard_name = "latitude_north"
        latf.axis = "y"
        latf.coordinate_defines = "center"
        latf[:]=np.float32(args[1])

    if cntrl >= 3:
        outf.createDimension('lon',args[2].shape[0])
        lonf = outf.createVariable('lon','f',('lon',))
        lonf.units = "degrees_east"
        lonf.long_name = "Longitude"
        lonf.actual_range = [args[2][0], args[2][-1]]
        #lonf.standard_name = "longitude_east"
        lonf.axis = "x"
        lonf.coordinate_defines = "center"
        lonf[:]=np.float32(args[2])

    if cntrl >= 4:
        outf.createDimension('level',args[3].shape[0])
        levf = outf.createVariable('level','f',('level',))
        levf.units = "millibar"
        levf.actual_range = [args[3][0], args[3][-1]]
        levf.long_name = "Level"
        levf.positive = "down"
        levf.GRIB_id = 100
        levf.GRIB_name = "hPa"
        levf.axis = "z"
        levf.coordinate_defines = "point"
        levf[:]=np.float32(args[3])

    if cntrl == 1:
        dataf = outf.createVariable(datastr,'f',('time'))
    elif cntrl == 2:
        dataf = outf.createVariable(datastr,'f',('time','lat'))
    elif cntrl == 3:
        dataf = outf.createVariable(datastr,'f',('time','lat','lon'))
    elif cntrl == 4:
        #SOMETIMES HAND EDIT IS NEEDED HERE DEPEND ON 
        #PERMUTATION OF YOUR DATA DIMENSIONS
        #dataf = outf.createVariable(datastr,'f',('time','lat','lon','level'))
        dataf = outf.createVariable(datastr,'f',('time','level','lat','lon'))
    dataf.missing_value = 32766.0
    #dataf.scale_factor = float(0.01)
    #dataf.long_name = "Sea Surface Temperature"
    #dataf.valid_range = [-5, 40]
    #dataf.actual_range = [np.where(np.isnan(data),0,data).min(), np.where(np.isnan(data),0,data).max()]
    #dataf.units = "degC"
    #dataf.add_offset = 0
    #dataf.precision = 2
    #dataf.least_significant_digit = 2
    #dataf.var_desc = "Sea Surface Temperature"
    #dataf.dataset = "NOAA Optimum Interpolation (OI) SST V2"
    #dataf.level_desc = "Surface"
    #dataf.statistic = "Mean"
    #dataf.parent_stat = "Weekly Mean"
    #dataf.standard_name = "sea_surface_temperature"
    #dataf.cell_methods = "time: mean (monthly from weekly values interpolated to daily)"
    try:
        dataf[:] = np.where(np.isnan(data),dataf.missing_value,data)
    except MemoryError:
        print "writenc: File too big, looping through time..."
        for it in xrange(data.shape[0]):
            dataf[it,:] = np.where(np.isnan(data[it,:]),\
                                   dataf.missing_value,data[it,:])
    #dataf[:]=data
    outf.close()

def writencep2(ncfile,datastr,False,data,*args,**kwargs):
    timeunits="hours since 1800-01-01 00:00"
    writenc(ncfile,datastr,timeunits,'ncep2',data,*args,**kwargs)

def writenc_generic(ncfile,datastr,data,*args,**kwargs):
    timeunits="hours since 1800-01-01 00:00"
    writenc(ncfile,datastr,timeunits,'ncep2',data,*args,**kwargs)

def openoisstclim(ncfile,varstr='sst'):
    '''Needed special function because this is long-term mean so no
       actual valid dates.'''
    dimlist = dimdict['ncep2'][:]
    dimlist.remove('level')
    dimlist.remove('lev')
    timestr = dimlist[0]
    ncf = kh.NetCDFFile(ncfile,'r')
    #vkeys = ncf.variables.keys() # could use something with this?
    for i in dimlist[:]:
            exec(i + ' = np.float32(ncf.variables[\''+ i + '\'][:])')
    
    exec('data = ncf.variables[\''+ varstr + '\'][:]')
    exec(varstr+' = addscalefill(ncf,data,varstr)')
    return sst, time, lat, lon
