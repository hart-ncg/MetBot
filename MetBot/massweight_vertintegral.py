import numpy as np

def massweight_vertintegral(var,levels,sfc_prs,prs_top,lev_units="hPa"):
    ''' integrated = massweight_vertintegral(var,levels,sfc_prs,prs_top,
                                            lev_units="hPa")
        Computed the mass weighted vertical integral of given gridded variable:
        verticalintegral = sum(var*dp/g) between surface and specified prs_top
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
    '''
    nax=np.newaxis
    g = 9.8
    nt,nz,ny,nx = var.shape
    ### Only loop in time for large arrays
    if var.nbytes/1024/1024 < 700:
        if (levels[0]-levels[1]) < 0:
            print "Levels need to be decreasing, inverting levels and variable"
            levels=levels[::-1]
            var = var[:,::-1,:,:]
        levsgrid = np.tile(levels[nax,:,nax,nax],(nt,1,ny,nx))
        sfcgrid = np.tile(sfc_prs[:,nax,:,:],(1,nz,1,1)) 
        topgrid = np.tile(np.float32(prs_top),(nt,nz,ny,nx))
        ilvs = np.int8(np.arange(len(levels)))
        ilvgrid = np.tile(ilvs[nax,:,nax,nax],(nt,1,ny,nx))
        validlevels = ~( ((levsgrid-topgrid)>0) & ((levsgrid-sfcgrid)<0) )
        levelsmasked = np.ma.MaskedArray(ilvgrid,mask = validlevels)
        ivalidbottom = levelsmasked.min(1).data
        ivalidbottom = np.tile(ivalidbottom[:,nax,:,:],(1,nz,1,1))
        lowestlev_mask = ivalidbottom==ilvgrid
        ivalidtop = levelsmasked.max(1).data
        ivalidtop = np.tile(ivalidtop[:,nax,:,:],(1,nz,1,1))
        upperlev_mask = ivalidtop==ilvgrid
        levs_c = np.float32((levels[1:] + levels[:-1])/2.)
        clevsgrid = np.tile(levs_c[nax,:,nax,nax],(nt,1,ny,nx))
        dp=np.zeros((nt,nz-1,ny,nx))
        dp[:,1:,:,:] = clevsgrid[:,:-1,:,:] - clevsgrid[:,1:,:,:]
        ### Now replace dp with with sfc dp in lowest valid level
        dplowest2rest=sfcgrid[:,:-1,:,:]-clevsgrid
        dphighest2rest=clevsgrid-topgrid[:,:-1,:,:]
        dp = np.where(lowestlev_mask[:,:-1,...],dplowest2rest[:,:,...],dp)
        dp = np.where(upperlev_mask[:,:-1,...],dphighest2rest[:,:,...],dp)
        massweighted_var = 100*var[:,:-1,:,:]*dp/g
        maskmwv=np.ma.MaskedArray(massweighted_var,mask=validlevels[:,:-1,...])
        vertical_integral = maskmwv.sum(1)
    else:
        print '''This is a big array, implement looping in time in this method
                 Doing nothing now.
              '''

    return vertical_integral.data

