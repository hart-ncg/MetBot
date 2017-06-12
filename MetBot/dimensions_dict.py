dim_deets={}

dim_deets['olr']={
"noaa": ['time','lat','lon'],
"ncep": ['time','lat','lon'],
"era": ['time','latitude','longitude'],
"20cr": ['time','lat','lon'],
"um": ['t','latitude','longitude','toa'],
"cmip5": ['time','lat','lon'],
}

dim_deets['pr']={
"trmm": ['time','latitude','longitude'],
"ncep": ['time','lat','lon'],
"era": ['time','latitude','longitude'],
"20cr": ['time','lat','lon'],
"um": ['t','latitude','longitude','surface'],
"cmip5": ['time','lat','lon'],
}

dim_deets['omega']={
"ncep": ['time','lat','lon','level'],
"era": ['time','latitude','longitude'],
"cmip5": ['time','lat','lon','plev'],
}

dim_deets['u']={
"ncep": ['time','lat','lon','level'],
"cmip5": ['time','lat','lon','plev'],
}

dim_deets['v']={
"ncep": ['time','lat','lon','level'],
"cmip5": ['time','lat','lon','plev'],
}