import numpy as np

def massweight_vertintegral(var,levels,sfc_prs,prs_top,lev_units="hPa"):
    ''' integrated = massweight_vertintegral(var,levels,sfc_prs,prs_top,
                                            lev_units="hPa")
        Computed the mass weighted vertical integral of given gridded variable:
        verticalintegral = sum(var*dp/g) between surface and specified prs_top
            #  RJ i think up to but not including prs_top?
        Implementation close to that of Grads vint, and in principle same as
        non-mass corrected integrals computed at CGD group, NCAR
        (http://www.cgd.ucar.edu/cas/catalog/newbudgets/index.html)
        Mass correction NOT implemented as yet.
    
    USAGE: var    - time x lev x lat x lon
           levels - levels coordinate
           sfc_prs- surface pressure - ?a field (time x lat x lon)?
           prs_top- user-specified top level to integrate to needs to
                    be at least lower (higher P) than topmost level
           lev_units[="hPa"] equivalently "mb", indicating units of lev coord.
    RETURNS: mass weighted vertical integral
    Written by NCG Hart July 2017
    Comments added by R James August 2017
    '''
    nax=np.newaxis # shortcut to add a new axis
    g = 9.8 # gravitational acceleration
    nt,nz,ny,nx = var.shape # dimensions of main variable
    ### Only loop in time for large arrays
#    if var.nbytes/1024/1024 < 700:
    if (levels[0]-levels[1]) < 0:
        print "Levels need to be decreasing, inverting levels and variable"
        levels=levels[::-1]
        var = var[:,::-1,:,:]
    levsgrid = np.tile(levels[nax,:,nax,nax],(nt,1,ny,nx)) # grid with dimensions of var
            #  at each level values are all set to the level itself i.e. 1000, 925, 850, 700, etc.
    sfcgrid = np.tile(sfc_prs[:,nax,:,:],(1,nz,1,1))  # copy the surface pressure input at all levels
    topgrid = np.tile(np.float32(prs_top),(nt,nz,ny,nx)) # grid with dimensions of var, but all values are prs_top e.g. 700
    ilvs = np.int8(np.arange(len(levels))) # integer counting each level, 0,1,2,3,4 etc. to nz in levels
    ilvgrid = np.tile(ilvs[nax,:,nax,nax],(nt,1,ny,nx)) # grid size of var but with these counts at each level
        # equivalent to levsgrid but with counts - 0s instead of 1000s, 1s instead of 925s etc.
    validlevels = ~( ((levsgrid-topgrid)>0) & ((levsgrid-sfcgrid)<0) ) # boolean grid True/False to show levels we want to integrate
            # removes levels above prs_top and below sfc_prs
            # but inverted so the ones we want are False (inverted using the tilde)
            # so it gives "True" for all the levels we want to mask
            # also it goes up to "prs_top" but not including
    levelsmasked = np.ma.MaskedArray(ilvgrid,mask = validlevels) # Masks levels we dont want to integrate (from "count" grid)
    ivalidbottom = levelsmasked.min(1).data # shape (time, lat, lon) for the count of the level we integrate from (often 0 or 1)
    ivalidbottom = np.tile(ivalidbottom[:,nax,:,:],(1,nz,1,1)) # copying this to all levels
    lowestlev_mask = ivalidbottom==ilvgrid # mask to show where masked because below surface
    ivalidtop = levelsmasked.max(1).data # shape (time, lat, lon) for the count of the level we integrate to
        # (will be count of level below prs_top)
    ivalidtop = np.tile(ivalidtop[:,nax,:,:],(1,nz,1,1)) # copying this to all levels
    upperlev_mask = ivalidtop==ilvgrid  # mask to show where masked because above prs_top
            # wait no this only shows "True" for the level below prs_top
    levs_c = np.float32((levels[1:] + levels[:-1])/2.) # hPa value for centre of each pressure box - 962.5, 887.5 etc.
    clevsgrid = np.tile(levs_c[nax,:,nax,nax],(nt,1,ny,nx)) # grid which copies these values across each level
        # similar to levsgrid but for centre of each pressure box
    dp=np.zeros((nt,nz-1,ny,nx)) # empty array with dimensions of var (except the '-1' 1 less for levels - excluding the top)
    dp[:,1:,:,:] = clevsgrid[:,:-1,:,:] - clevsgrid[:,1:,:,:] # difference in pressure between each level
        # the 1: misses the bottom level
    ### Now replace dp with with sfc dp in lowest valid level
    dplowest2rest=sfcgrid[:,:-1,:,:]-clevsgrid
    dphighest2rest=clevsgrid-topgrid[:,:-1,:,:] # diff to prs_top
    dp = np.where(lowestlev_mask[:,:-1,...],dplowest2rest[:,:,...],dp)
        # where there is the lowest lev mask replace with the value which calcs diff to surface
    dp = np.where(upperlev_mask[:,:-1,...],dphighest2rest[:,:,...],dp)
        # where there is the upper lev mask (only on level below prs_top) replace with the value which calcs diff to surface
    massweighted_var = 100*var[:,:-1,:,:]*dp/g # calculates mass weighted var for each gridpoint
        # multiplication by 100 explains why the values are larger than I expected
    maskmwv=np.ma.MaskedArray(massweighted_var,mask=validlevels[:,:-1,...]) # masks values specified above
    vertical_integral = maskmwv.sum(1) # sum over remaining mass weighted values
#    else:
#        print '''This is a big array, implement looping in time in this method
#                 Doing nothing now.
#              '''

    return vertical_integral.data


def massweight_vint_nosurf(var, levels, prs_top, lev_units="hPa"):
    ''' integrated = massweight_vertintegral(var,levels,prs_top,
                                            lev_units="hPa")
        Computed the mass weighted vertical integral of given gridded variable:
        verticalintegral = sum(var*dp/g) between surface and specified prs_top
            #  RJ i think up to but not including prs_top?
        Works on files where the values below surface are already set to missing
        e.g. some CMIP5 models
        Implementation close to that of Grads vint, and in principle same as
        non-mass corrected integrals computed at CGD group, NCAR
        (http://www.cgd.ucar.edu/cas/catalog/newbudgets/index.html)
        Mass correction NOT implemented as yet.

    USAGE: var    - time x lev x lat x lon
           levels - levels coordinate
           sfc_prs- surface pressure - ?a field (time x lat x lon)?
           prs_top- user-specified top level to integrate to needs to
                    be at least lower (higher P) than topmost level
           lev_units[="hPa"] equivalently "mb", indicating units of lev coord.
    RETURNS: mass weighted vertical integral
    Written by NCG Hart July 2017
    Comments added by R James August 2017
    '''
    nax = np.newaxis  # shortcut to add a new axis
    g = 9.8  # gravitational acceleration
    nt, nz, ny, nx = var.shape  # dimensions of main variable
    ### Only loop in time for large arrays
    #    if var.nbytes/1024/1024 < 700:
    if (levels[0] - levels[1]) < 0:
        print "Levels need to be decreasing, inverting levels and variable"
        levels = levels[::-1]
        var = var[:, ::-1, :, :]
    levsgrid = np.tile(levels[nax, :, nax, nax], (nt, 1, ny, nx))  # grid with dimensions of var
    #  at each level values are all set to the level itself i.e. 1000, 925, 850, 700, etc.
    topgrid = np.tile(np.float32(prs_top),
                      (nt, nz, ny, nx))  # grid with dimensions of var, but all values are prs_top e.g. 700
    ilvs = np.int8(np.arange(len(levels)))  # integer counting each level, 0,1,2,3,4 etc. to nz in levels
    ilvgrid = np.tile(ilvs[nax, :, nax, nax], (nt, 1, ny, nx))  # grid size of var but with these counts at each level
    # equivalent to levsgrid but with counts - 0s instead of 1000s, 1s instead of 925s etc.
    validlevels = ~((levsgrid - topgrid) > 0) # boolean grid True/False to show levels we want to integrate
    # removes levels above prs_top
    # but inverted so the ones we want are False (inverted using the tilde)
    # so it gives "True" for all the levels we want to mask
    # also it goes up to "prs_top" but not including
    levelsmasked = np.ma.MaskedArray(ilvgrid,
                                     mask=validlevels)  # Masks levels we dont want to integrate (from "count" grid)
    ivalidbottom = levelsmasked.min(1).data  # shape (time, lat, lon) for the count of the level we integrate from (often 0 or 1)
    ivalidbottom = np.tile(ivalidbottom[:, nax, :, :], (1, nz, 1, 1))  # copying this to all levels
    lowestlev_mask = ivalidbottom == ilvgrid  # mask to show where masked because below surface
    ivalidtop = levelsmasked.max(1).data  # shape (time, lat, lon) for the count of the level we integrate to
    # (will be count of level below prs_top)
    ivalidtop = np.tile(ivalidtop[:, nax, :, :], (1, nz, 1, 1))  # copying this to all levels
    upperlev_mask = ivalidtop == ilvgrid  # mask to show where masked because above prs_top
    # wait no this only shows "True" for the level below prs_top
    levs_c = np.float32(
        (levels[1:] + levels[:-1]) / 2.)  # hPa value for centre of each pressure box - 962.5, 887.5 etc.
    clevsgrid = np.tile(levs_c[nax, :, nax, nax], (nt, 1, ny, nx))  # grid which copies these values across each level
    # similar to levsgrid but for centre of each pressure box
    dp = np.zeros((nt, nz - 1, ny,
                   nx))  # empty array with dimensions of var (except the '-1' 1 less for levels - excluding the top)
    dp[:, 1:, :, :] = clevsgrid[:, :-1, :, :] - clevsgrid[:, 1:, :, :]  # difference in pressure between each level
    # the 1: misses the bottom level
    ### Now replace dp with with sfc dp in lowest valid level
    dplowest2rest = sfcgrid[:, :-1, :, :] - clevsgrid
    dphighest2rest = clevsgrid - topgrid[:, :-1, :, :]  # diff to prs_top
    dp = np.where(lowestlev_mask[:, :-1, ...], dplowest2rest[:, :, ...], dp)
    # where there is the lowest lev mask replace with the value which calcs diff to surface
    dp = np.where(upperlev_mask[:, :-1, ...], dphighest2rest[:, :, ...], dp)
    # where there is the upper lev mask (only on level below prs_top) replace with the value which calcs diff to surface
    massweighted_var = 100 * var[:, :-1, :, :] * dp / g  # calculates mass weighted var for each gridpoint
    # multiplication by 100 explains why the values are larger than I expected
    maskmwv = np.ma.MaskedArray(massweighted_var, mask=validlevels[:, :-1, ...])  # masks values specified above
    vertical_integral = maskmwv.sum(1)  # sum over remaining mass weighted values
    #    else:
    #        print '''This is a big array, implement looping in time in this method
    #                 Doing nothing now.
    #              '''

    return vertical_integral.data
