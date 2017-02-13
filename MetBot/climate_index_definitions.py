#!/usr/bin/env python
# This contains a dictionary of lat x lon definitions of well know climate indices
# Dict entries follow format index_key: [latmin,latmax,lonmin,lonmax]
# lats convection -90:90, lons 0:359

indexdict={}
indexdict['nino12'] = [-10,0,270,280]
indexdict['nino3']  = [-5,5,210,270]
indexdict['nino34'] = [-5,5,190,240]
indexdict['nino4']  = [-5,5,160,210]
