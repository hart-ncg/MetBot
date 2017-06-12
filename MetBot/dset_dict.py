# DICTIONARY FOR DATASETS
# Funtionality for multiple datasets including noaa, um, cmip5
# and now ncep, era, 20cr
# Structure $dset $model
# start date is the first date in the file as saved

dset_deets={}

# NOAA
dset_deets['noaa']={
'noaa': {'name': "noaa", 'calendar':"gregorian", 'olrtimeunit':"hours since 1800-01-01 00:00:0.0", 'prtimeunit':"hours since 1800-01-01 00:00:0.0", 'startdate':"1974-06-01", 'olrname':"olr", 'prname':"none", 'yrfname':"1974_2013", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
'cdr': {'name': "cdr", 'calendar':"gregorian", 'olrtimeunit':"days since 1800-01-01 00:00:0.0", 'prtimeunit':"days since 1800-01-01 00:00:0.0", 'startdate':"1979-01-01", 'olrname':"olr", 'prname':"none", 'yrfname':"1979_2012", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
}

# TRMM
dset_deets['trmm']={
'trmm_3b42v7': {'name': "trmm_3b42v7", 'calendar':"standard", 'olrtimeunit':"hours since 1800-01-01 00:00:0.0", 'prtimeunit':"hours since 1998-1-1 03:00:0.0", 'startdate':"1998-01-01", 'olrname':"none", 'prname':"pcp", 'yrfname':"1998_2013", 'startyr':"1998", 'testfileyr':"1998_1998",'testyr':"1998"},
}

# NCEP
dset_deets['ncep']={
'ncep2': {'name': "ncep2", 'calendar':"gregorian", 'olrtimeunit':"hours since 1800-01-01 00:00:0.0", 'prtimeunit':"hours since 1800-01-01 00:00:0.0", 'omegatimeunit':"hours since 1800-1-1 00:00:00", 'startdate':"1979-01-01", 'olrname':"ulwrf", 'prname':"prate", 'yrfname':"1979_2013", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
}

# ERA
dset_deets['era']={
'erai': {'name': "erai", 'calendar':"gregorian", 'olrtimeunit':"hours since 1900-01-01 00:00:0.0", 'prtimeunit':"hours since 1900-01-01 00:00:0.0", 'startdate':"1979-01-01", 'olrname':"olr", 'prname':"tp", 'yrfname':"1979_2013", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
}

# 20cr
dset_deets['20cr']={
'20CRv2c': {'name': "20CRv2c", 'calendar':"standard", 'olrtimeunit':"hours since 1800-01-01 00:00:0.0", 'prtimeunit':"hours since 1800-01-01 00:00:0.0", 'startdate':"1978-01-01", 'olrname':"ulwrf", 'prname':"prate", 'yrfname':"1978_2012", 'startyr':"1978", 'testfileyr':"1979_1979",'testyr':"1979"},
}

# UM
dset_deets['um']={
'anqjn': {'name': "anqjn", 'calendar':"360_day", 'olrtimeunit':"days since 1978-09-01 00:00:00", 'prtimeunit':"days since 1978-09-01 00:00:00", 'startdate':"1978-12-01", 'olrname':"olr", 'prname':"precip", 'yrfname':"1978_2013", 'startyr':"1979", 'testfileyr':"1979_1979",'testyr':"1979"},
'antib': {'name': "antib", 'calendar':"360_day", 'olrtimeunit':"days since 1981-09-01 00:00:00", 'prtimeunit':"days since 1981-09-01 00:00:00", 'startdate':"1981-09-01", 'olrname':"olr", 'prname':"precip", 'yrfname':"1981_2008", 'startyr':"1982", 'testfileyr':"1985_1985",'testyr':"1985"},
'u-ab674': {'name': "u-ab674", 'calendar':"360_day", 'olrtimeunit':"days since 1998-12-1 00:00:00", 'prtimeunit':"days since 1998-12-1 00:00:00", 'startdate':"2079-01-01", 'olrname':"olr", 'prname':"precip", 'yrfname':"2079_2113", 'startyr':"2079", 'testfileyr':"2079_2079",'testyr':"2079"},
'u-ab680': {'name': "u-ab680", 'calendar':"360_day", 'olrtimeunit':"days since 1981-9-1 00:00:00", 'prtimeunit':"days since 1981-9-1 00:00:00", 'startdate':"1981-09-01", 'olrname':"olr", 'prname':"precip", 'yrfname':"1981_2008", 'startyr':"1982", 'testfileyr':"1985_1985",'testyr':"1985"},
}

# CMIP5
dset_deets['cmip5']={
'ACCESS1-0': {'name': "ACCESS1-0", 'calendar':"proleptic_gregorian", 'olrtimeunit':"days since 1-01-01 00:00:00", 'prtimeunit':"days since 1-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'ACCESS1-3': {'name': "ACCESS1-3", 'calendar':"proleptic_gregorian", 'olrtimeunit':"days since 1-01-01 00:00:00", 'prtimeunit':"days since 1-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'bcc-csm1-1-m': {'name': "bcc-csm1-1-m", 'calendar':"365_day", 'olrtimeunit':"days since 1950-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'BNU-ESM': {'name': "BNU-ESM", 'calendar':"365_day", 'olrtimeunit':"days since 1950-01-01 00:00:00", 'prtimeunit':"days since 1950-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CanESM2': {'name': "CanESM2", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CMCC-CESM': {'name': "CMCC-CESM", 'calendar':"standard", 'olrtimeunit':"days since 1950-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00",'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CMCC-CM': {'name': "CMCC-CM", 'calendar':"standard", 'olrtimeunit':"days since 1950-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CMCC-CMS': {'name': "CMCC-CMS", 'calendar':"standard", 'olrtimeunit':"days since 1950-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CNRM-CM5': {'name': "CNRM-CM5", 'calendar':"standard", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'CSIRO-Mk3-6-0': {'name': "CSIRO-Mk3-6-0", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'FGOALS-g2': {'name': "FGOALS-g2", 'calendar':"365_day", 'olrtimeunit':"days since 1-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'GFDL-CM3': {'name': "GFDL-CM3", 'calendar':"365_day", 'olrtimeunit':"days since 1860-01-01 00:00:00", 'prtimeunit':"days since 1860-01-01 00:00:00", 'omegatimeunit':"days since 1860-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'omeganame':"wap", 'fullrun':"1860_2005", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'GFDL-ESM2G': {'name': "GFDL-ESM2G", 'calendar':"365_day", 'olrtimeunit':"days since 1861-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'GFDL-ESM2M': {'name': "GFDL-ESM2M", 'calendar':"365_day", 'olrtimeunit':"days since 1861-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'HadGEM2-CC': {'name': "HadGEM2-CC", 'calendar':"360_day", 'olrtimeunit':"days since 1859-12-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'HadGEM2-ES': {'name': "HadGEM2-ES", 'calendar':"360_day", 'olrtimeunit':"days since 1859-12-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'inmcm4': {'name': "inmcm4", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'IPSL-CM5A-LR': {'name': "IPSL-CM5A-LR", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'IPSL-CM5A-MR': {'name': "IPSL-CM5A-MR", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'IPSL-CM5B-LR': {'name': "IPSL-CM5B-LR", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MIROC5': {'name': "MIROC5", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MIROC-ESM-CHEM': {'name': "MIROC-ESM-CHEM", 'calendar':"standard", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MIROC-ESM': {'name': "MIROC-ESM", 'calendar':"standard", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'omegatimeunit':"days since 1850-1-1 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'omeganame':"wap", 'fullrun':"1950_2005", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MPI-ESM-LR': {'name': "MPI-ESM-LR", 'calendar':"proleptic_gregorian", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MPI-ESM-MR': {'name': "MPI-ESM-MR", 'calendar':"proleptic_gregorian", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'MRI-CGCM3': {'name': "MRI-CGCM3", 'calendar':"standard", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
'NorESM1-M': {'name': "NorESM1-M", 'calendar':"365_day", 'olrtimeunit':"days since 1850-01-01 00:00:00", 'prtimeunit':"days since 1850-01-01 00:00:00", 'startdate':"1970-01-01", 'olrname':"rlut", 'prname':"pr", 'yrfname':"1970_2004", 'startyr':"1970", 'testfileyr':"1975_1975",'testyr':"1975"},
}