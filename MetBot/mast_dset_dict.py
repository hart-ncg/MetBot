# DICTIONARY FOR DATASETS (not models)
# this is for things that are consistent across datasets

mast_dset_deets={}

# NOAA
mast_dset_deets['noaa']={
'olrname':"olr",
}

# TRMM
dset_deets['trmm']={
'prname':"pcp",
}

# NCEP
mast_dset_deets['ncep']={
'olrname':"ulwrf",'prname':"prate",'omeganame':"omega",'uname':"uwnd",'vname':"vwnd", 'qname':"shum", 'gpthname':"hgt", 'Tname':"air", 'presname':"pres",
}

# ERA
dset_deets['era']={
'olrname':"olr", 'prname':"tp", 'omeganame':"w", 'uname':"u", 'vname':"v", 'qname':"q", 'gpthname':"z", 'Tname':"t", 'presname':"sp",
}

# CMIP5
mast_dset_deets['cmip5']={
'olrname':"rlut", 'prname':"pr", 'omeganame':"wap", 'uname':"ua", 'vname':"va", 'qname':"hus", 'gpthname':"zg", 'Tname':"ta",
}

