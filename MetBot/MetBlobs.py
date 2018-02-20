# MetBlobs.py
''' MetBlobs.py Module: Rudimentary blob detection based feature identification system
cvblob class has attributes: minx, maxx, miny, centroid(x,y), area,
a metblobs record (mbs) has first dim time and then 14 values as below
HOURtime, Label, Angle, cX, cY, minX, maxX, minY, maxY, area, Circularity, varmean_in_convhull, vrbmn_in_convhull, vrbmx_in_convhull
'''
import Image
import skimage, skimage.color, skimage.measure, skimage.io, skimage.draw
import skimage.morphology
import os
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.nxutils as nx ### Obselete since matplotlib 1.3
from matplotlib.path import Path ### This implements updated nx
import mpl_toolkits.basemap as bm
try: from mpl_toolkits.basemap import num2date
except ImportError: from netCDF4 import num2date
#import mytools as my
from time import time as timer
import cPickle
debugfilters=False
# SOME DEFAULT PARAMETERS
# DICTIONARY OF FILTER VALUES
import MetBot.filterdict_default as filters

# THRESHOLD LIMITS SET UP
def calclims(vrb,scalestd=1.3,disttype='normal'):
    '''Calculates threshold appropriate to given distributions
    Options for type: normal or bimodal
    for normal dist, stretch vals are those at the 99th percentile (two-tailed)
    thresh's are specified as a value of the standard deviation'''
    if disttype=='normal':
        vm = vrb.mean()
        vs = vrb.std()
        #xg=sp.stats.norm.pdf(x, loc=vm, scale=vs)
        tlow, thigh = vm-scalestd*vs, vm+scalestd*vs
        slow, shigh = vm-3*vs, vm+3*vs
    return (slow, shigh), (tlow, thigh)

# HELPER FUNCTIONS
def mbsave(fname,mbs,mbt,chull):
    print "Saving metblobs to", fname
    pickf = open(fname,'w')
    cPickle.dump((mbs,mbt,chull),pickf)
    pickf.close()

def mbopen(fname):
    '''mbs, mbt, chull = mbopen(fname)'''
    print "Opening", fname
    pickf = open(fname,'r')
    mbs, mbt, chull = cPickle.load(pickf)
    pickf.close()
    return mbs, mbt, chull

def mbtwrite(fname,mbt,mbs):
    '''Simple metblobs time to text file writer'''
    print "Writing metblobs time to", fname
    pickf = open(fname,'w')
    for i in xrange(len(mbt)):
        t=mbt[i]
        pickf.write("%d,%02d,%02d,%02d \n" %(t[0],t[1],t[2],np.int(mbs[i,3])))
    pickf.close()

def mmm_ch_val(vrb,llgridtup,ch):
    '''Function to calculate mean of met variable inside chull'''
    ln,lt = llgridtup
    lls = ln.shape
    xypts = np.hstack((ln.ravel()[:,np.newaxis], lt.ravel()[:,np.newaxis]))
    tes = np.reshape(xypts[:,0],lls)
    chpoly = Path(ch)
    mask = chpoly.contains_points(xypts)
    maskgrid = np.reshape(mask, lls)
    maskgrid = np.where(maskgrid,False,True)
    vrbma=np.ma.MaskedArray(vrb,mask=maskgrid)

    return vrbma.mean(), vrbma.min(), vrbma.max()

def randindgen(idxlen,n_irnd):
    '''irandom = randindgen(idxlen,n_irnd)
    Random Index Generator
    Generates a set of random integers on the interval [0,idxlen)
    Useful as index of random slices from a large array'''
    irandom = np.int16(np.round(np.random.random((n_irnd))*(idxlen-1)))
    return  irandom

class Geopix:
    '''Geopix((npix_x,npix_y),lat,lon)
    Convert between image pixels and geolocations
    Methods:
    lat2pix(deg)
    lat2pix(deg)
    imgsubdomainmasked(img,domain)
    xp2lon(xpix)
    yp2lat(ypix)'''

    def __init__(self,(npix_y,npix_x),lat,lon):
        self.rnglon = lon[-1]-lon[0]
        self.rnglat = lat[-1]-lat[0]
        self.lon0=lon[0]
        self.lat0=lat[0]
        self.lat=lat;self.lon=lon
        self.nx=npix_x;self.ny=npix_y

    def lat2pix(self,deg):
        #ideg=np.where(np.abs(self.lat-deg)==np.abs((self.lat-deg)).min())[0][0]
        #closedeg = self.lat[ideg]
        pix=(deg-self.lat0)*self.ny/float(self.rnglat)
        #print pix
        return np.int(pix)

    def lon2pix(self,deg):
        #ideg=np.where(np.abs(self.lon-deg)==np.abs((self.lon-deg)).min())[0][0]
        #closedeg = self.lon[ideg]
        pix=(deg-self.lon0)*self.nx/float(self.rnglon)
        #print 'lon2pix',pix,self.lon[ideg]
        return np.int(pix)

    #def imgsubdomain(self,img,domain):
    #    if not domain:
    #        lon1,lon2 = self.lon[0], self.lon[-1]
    #        lat1,lat2 = self.lat[0], self.lat[-1]
    #    else:
    #        lon1,lon2,lat1,lat2 = domain
    #    x1,x2 = self.lon2pix(lon1),self.lon2pix(lon2)
    #    y1,y2 =  self.lat2pix(lat1),self.lat2pix(lat2)
    #    subrect=(y1,x1,y2-y1,x2-x1)
    #    x1,y1,w,h=subrect
    #    if img.ndim>2: subimg = img[x1:x1+w,y1:y1+h,:]
    #    elif img.ndim==2: subimg = img[x1:x1+w,y1:y1+h]
    #    self.subdomain=(x1,x2,y1,y2,w,h)

    #    return subimg, (x1,x2,y1,y2,w,h)

    def imgsubdomainmasked(self,img,domain):
        if not domain:
            lon1,lon2 = self.lon[0], self.lon[-1]
            lat1,lat2 = self.lat[0], self.lat[-1]
        else:
            lon1,lon2,lat1,lat2 = domain
        x1,x2 = self.lon2pix(lon1),self.lon2pix(lon2)
        y1,y2 =  self.lat2pix(lat1),self.lat2pix(lat2)
        subrect=(y1,x1,y2-y1,x2-x1)
        #x1,y1,w,h=subrect
        y1,x1,w,h=subrect
        #print img.shape
        yy,xx = np.indices(img.shape)
        submask = ((xx>x1) & (xx<x2)) & ((yy>y1) & (yy<y2))
        black=0
        if img.ndim>2:subimg = np.where(submask,img,black)
        elif img.ndim==2: subimg = np.where(submask,img,0)
        self.subdomain=(x1,x2,y1,y2,w,h)

        return subimg, (x1,x2,y1,y2,w,h)

    def xp2lon(self,xpix):
        #x1 = self.subdomain[0]
        #pix = xpix+x1
        pix = xpix
        ln = pix*self.rnglon/float(self.nx) + self.lon0
        return np.around(ln,decimals=1)

    def yp2lat(self,ypix):
        #y1 = self.subdomain[2]
        #pix = ypix+y1
        pix = ypix
        lt = pix*self.rnglat/float(self.ny) + self.lat0
        return np.around(lt,decimals=1)

#def fig2img(fig,subrect):
#    '''
#    Convert a Matplotlib figure to a ndarray image in RGB format and return it
#    @param fig a matplotlib figure
#    @return a ndarray image
#    Adapted from: http://www.icare.univ-lille1.fr/wiki/index.php/How_to_convert_a_matplotlib_figure_to_a_numpy_array_or_a_PIL_image'''
#    # draw the renderer
#    fig.canvas.draw()
#    # Get the RGBA buffer from the figure
#    w,h = fig.canvas.get_width_height()
#    buf = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
#    img = buf.reshape((h,w,3))
#    x1,y1,xw,yh=subrect
#    if img.ndim>2: plotarea = img[x1:x1+xw,y1:y1+yh,:]
#    elif img.ndim==2: plotarea = img[x1:x1+xw,y1:y1+yh]
#
#    return plotarea

#def CanvasCutOut(fig):
#    ''' subrect = CanvasCutout(fig)
#    A function that identifies the plotarea in a figure canvas and returns
#    the sub-rectangle that identifies this plotarea region.'''
#    fig.set_facecolor('k')
#    fig.canvas.draw()
#    w,h = fig.canvas.get_width_height()
#    buf = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
#    ndimg = buf.reshape((h,w,3))
#    ndimg_g = skimage.img_as_ubyte(skimage.color.rgb2gray(ndimg))
#    labelled = skimage.measure.label(ndimg_g>100,connectivity=2)
#    blobs = skimage.measure.regionprops(labelled)
#    sizes = np.asarray([b.area for b in blobs])
#    bigblob = blobs[sizes.argmax()]
#    x1, y1, x2, y2 = np.int32(bigblob.bbox)
#    xw, yh = x2-x1, y2-y1
#    subrect=(x1,y1,xw,yh)
#
#    return subrect

#def greatestblob(blobs):
#    '''Return the greater blob and its index in the blob dict'''
#    bk=blobs.keys()
#    bk.sort()
#    sizes = np.asarray([blobs[b].area for b in bk])
#    ixbig = sizes.argmax()
#    bigblob = blobs[ixbig]
#    return bigblob, ixbig

#def Stretch(arr,dthresh,fn="Linear",tval=None):
#    '''stretched = Stretch(arr,fn="Linear",tval=None)
#    tval = (threshmin,thresmax) stretches between user defined limits
#    default tval=None stretches between min&max of array'''
#    if tval:
#        #print "Stretching data to 0 255..."
#        arrIn=arr.copy()
#        arrOut = np.where(arrIn < tval[0],tval[0],arrIn)
#        arrOut = np.where(arrIn > tval[1],tval[1],arrOut)
#    arrOut=arrOut-tval[0]
#    m = 255.0/(tval[1]-tval[0])
#    greythresh=(dthresh-tval[0])*m
#    return np.round(arrOut*m), np.round(greythresh)
#
#def Threshold(arr,tval):
#    '''Black-White thresholding'''
#    arrIn=arr.copy()
#    arrOut = np.where(arrIn < tval,0,arrIn)
#    arrOut = np.where(arrIn > tval,255,arrOut)
#    return arrOut


def MagDir(u,v):
    '''returns the mag and bearing of vector describe by u & v'''
    mag=np.sqrt(u**2 + v**2)
    raddir=np.arctan(v/u)
    degdir=450-360*(raddir/(2*np.pi))
    return mag, raddir

def DistVector(p1,p2):
    ''' (x,y), (mag, dir) = DistVector(p1,p2)
    Input: p1 = np.array(x1, y1) & p2 = np.array(x2, y2)
    Returns the distance vector in forms
    (x,y) tuple
    (mag,dir) tuple'''
    x1=p1[:,0];y1=p1[:,1]
    x2=p2[:,0];y2=p2[:,1]
    dx=x2-x1
    dy=y2-y1
    mag,dir=MagDir(dx,dy)

    return (dx,dy),(mag,dir)

#def BlobAngles(dct,blobs,img,gpx,thrs):
def BlobAngles(dct,blobs,img,gpx):
    if 'angleROI' in dct.keys():   # Use an angle from blobs in subdomain of one used to get blobs passed to FilterByAngle
        angledomain = dct['angleROI']
        #bigsubdomain = gpx.subdomain
        img4blob, pixss = gpx.imgsubdomainmasked(img,angledomain)
        #smallsubdomain = gpx.subdomain
        #angblbim, angblobs, anggrey = GetBlobs(img4blob,(thrs,255))
        angblbim, angblobs, anggrey = GetBlobs(img4blob)
        #if plt.fignum_exists(3):
        #    lat, lon = gpx.lat, gpx.lon
        #    dlat=lat[1]-lat[0];dlon=lon[1]-lon[0]
        #    latplot = np.hstack((lat[0]-dlat/2.,lat+dlat/2.))
        #    lonplot = np.hstack((lon[0]-dlon/2.,lon+dlon/2.))
        #    plt.figure(num=3)
        #    plt.pcolormesh(lonplot,latplot,angblbim);plt.grid()
        #    plt.xlim(lonplot[0],lonplot[-1]);plt.ylim(latplot[-1],latplot[0])
        #    plt.draw()
        dstold = 100000 # Big irrelevant starting dst
        for i in blobs.keys():
            #gpx.subdomain=bigsubdomain
            xb = gpx.xp2lon(blobs[i].centroid[1])
            yb = gpx.yp2lat(blobs[i].centroid[0])
            ij = 1
            ### Associate angleblobs with blobs
            #print angblobs.keys()
            for j in angblobs.keys():
                #gpx.subdomain=smallsubdomain
                xab = gpx.xp2lon(angblobs[j].centroid[1])
                yab = gpx.yp2lat(angblobs[j].centroid[0])
                # should be distance on a sphere but isn't
                dst = np.sqrt((xab-xb)**2 + (yab-yb)**2)
                ### Test to get angleblob centroid closest to blob centroid
                if dst < dstold: ijmatch = j
            degs=(blobs[i].orientation*360)/(2*np.pi)
            angdegs=(angblobs[ijmatch].orientation*360)/(2*np.pi)
            blobs[i].degs = degs
            blobs[i].angdegs = angdegs
        #gpx.subdomain=bigsubdomain
    else:          # Nice and simple unlike above
        for i in blobs.keys():
            degs=(blobs[i].orientation*360)/(2*np.pi)
            blobs[i].degs = degs

def FilterByAngle(dct,blobs):
    minAngle, maxAngle = dct['angle']
    if 'angleROI' in dct.keys(): 
        for i in blobs.keys():
            ### Delete parent blob if angleblob angle fails
            if ( blobs[i].degs<minAngle or blobs[i].degs>maxAngle) & \
               ( blobs[i].angdegs<minAngle or blobs[i].angdegs>maxAngle): 
                del blobs[i]
                print "Failed angles"
    else:
        for i in blobs.keys():
            if ( blobs[i].degs<minAngle or blobs[i].degs>maxAngle):
                del blobs[i]

def FilterByLatextent(dct,blobs,gpx):
    latN,latS=dct['latextent']
    pixN, pixS = gpx.lat2pix(latN), gpx.lat2pix(latS)
    if hasattr(gpx,'subdomain'):
        y1, y2 = gpx.subdomain[0], gpx.subdomain[1]
        if pixN < y1:
            pixN = y1;
        if pixS > y2:
            pixS = y2;

    #mn, mx = pixN+5, pixS-5
    mn, mx = pixN, pixS
    #print "Lat pixs:", pixN, pixS
    for i in blobs.keys():
        #print 'Current iteration',i, blobs[i].label
        miny,minx,maxy,maxx = blobs[i].bbox
        #print miny, maxy, mn,mx
        if (miny>mn or maxy< mx):
            #print 'Deleting blob no.',i,' labelled:',blobs[i].label
            del blobs[i]
            #print "Failed Latextent"

#def FilterBlobs(dct,blobs,img,gpx,thrs):
def FilterBlobs(dct,blobs,img,gpx):
    '''Filters blobs by predefined criteria:
    These are contained in filters.blobfilters (dict type)'''
    # Area Filter
    if dct['area']:
        mn,mx=dct['area']
        #mn, mx = 4,300
        for i in blobs.keys():
            if blobs[i].area<mn or blobs[i].area>mx:
                #print 'Area:',blobs[i].area, "fails"
                del blobs[i]
    # Latitutdinal Extent
    if dct['latextent']:
        FilterByLatextent(dct,blobs,gpx)
    # Angle Filter
    #BlobAngles(dct,blobs,img,gpx,thrs)
    BlobAngles(dct,blobs,img,gpx)
    if dct['angle']:
        FilterByAngle(dct,blobs)



# MAPPING AND PLOTTING HELPERS

def SAfrBasemap(lat,lon,drawstuff=False,prj='cyl',fno=1,rsltn='c',\
    fontdict=False):
    '''m, f = SAfrBasemap(lat,lon,drawstuff=False,prj='cyl',fno=1)

    This creates a basemap instance from lat's and lon's provided.
    Specific for this application of metblobs, so currently using proj='cyl',
    however Albers equal area is here uncommented.
    USAGE: lat, lon
    RETURNS: m, basemap object pointer (handle in matlab language)
             f, pointer figure'''
    xy=(lat.min(),lat.max(),lon.min(),lon.max())
    nx = len(lon); ny = len(lat)
    if not fontdict: fontdict = {'fontsize':14,'fontweight':'bold'}
    if prj=='cyl':
        m = bm.Basemap(llcrnrlon=xy[2],llcrnrlat=xy[0],urcrnrlon=xy[3],\
                       urcrnrlat=xy[1],resolution=rsltn,area_thresh=10000.,\
                       projection='cyl')
    if prj=='aea':
        m = bm.Basemap(llcrnrlon=xy[2]-5.,llcrnrlat=xy[0],urcrnrlon=xy[3],\
                       urcrnrlat=xy[1]+5,resolution=rsltn,projection='aea',\
                        lat_1=-45.,lat_2=0.,lon_0=40.)

    ### SET UP FIGURE AND MERIDIANS
    delon = 10.
    meridians = np.arange(10.,360.,delon)
    delat = 5.
    circles = np.arange(0.,90.+delat,delat).tolist()+\
              np.arange(-delat,-90.-delat,-delat).tolist()
    f1 = plt.figure(fno,figsize=[13.0, 8.0])
    #ax = f1.add_axes([0.1,0.1,0.7,0.7])
    #pos = ax.get_position()
    #l, b, w, h = pos.bounds
    #cax = plt.axes([l+w+0.075, b, 0.05, h]) # setup colorbar axes.
    #plt.colorbar(cax=cax) # draw colorbar
    #plt.axes(ax)

    if drawstuff:
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(circles,linewidth='0.1',labels=[1,0,0,0],\
                        fontdict=fontdict)
        m.drawmeridians(meridians,linewidth='0.1',labels=[0,0,0,1],\
                        fontdict=fontdict)

    return m, f1

def SAfrBasemap2(lat,lon,drawstuff=False,prj='cyl', rsltn='c',\
    fontdict=False):
    '''m, f = SAfrBasemap(lat,lon,drawstuff=False,prj='cyl',fno=1)

    This creates a basemap instance from lat's and lon's provided.
    Specific for this application of metblobs, so currently using proj='cyl',
    however Albers equal area is here uncommented.
    USAGE: lat, lon
    RETURNS: m, basemap object pointer (handle in matlab language)
             f, pointer figure'''
    xy=(lat.min(),lat.max(),lon.min(),lon.max())
    nx = len(lon); ny = len(lat)
    if not fontdict: fontdict = {'fontsize':14,'fontweight':'bold'}
    if prj=='cyl':
        m = bm.Basemap(llcrnrlon=xy[2],llcrnrlat=xy[0],urcrnrlon=xy[3],\
                       urcrnrlat=xy[1],resolution=rsltn,area_thresh=10000.,\
                       projection='cyl')
    if prj=='aea':
        m = bm.Basemap(llcrnrlon=xy[2]-5.,llcrnrlat=xy[0],urcrnrlon=xy[3],\
                       urcrnrlat=xy[1]+5,resolution=rsltn,projection='aea',\
                        lat_1=-45.,lat_2=0.,lon_0=40.)

    ### SET UP FIGURE AND MERIDIANS
    delon = 10.
    meridians = np.arange(10.,360.,delon)
    delat = 5.
    circles = np.arange(0.,90.+delat,delat).tolist()+\
              np.arange(-delat,-90.-delat,-delat).tolist()

    if drawstuff:
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(circles,linewidth='0.1',labels=[1,0,0,0],\
                        fontdict=fontdict)
        m.drawmeridians(meridians,linewidth='0.1',labels=[0,0,0,1],\
                        fontdict=fontdict)

    return m

def BlobContour(blbim,blob,gpx):
    '''chtype=convexhull or simple (which is more complex)'''
    boolblob = blbim==blob.label
    polylist = skimage.measure.find_contours(boolblob,.5,\
               fully_connected='high',positive_orientation='high')
    if len(polylist)>1:
        #print "More than one poly Damn!"
        iply = np.asarray([len(p) for p in polylist]).argmax()
        poly=polylist[iply]
        #poly = np.vstack((poly,poly[0,:]))
        apoly=skimage.measure.approximate_polygon(poly,1)
        #nanbuf = np.asarray((np.nan,np.nan))[np.newaxis,:]
        #for ply in polylist[1:]:
        #    lastpnt=ply[0,:][np.newaxis,:]
        #    aply=skimage.measure.approximate_polygon(ply,1)
        #    apoly = np.vstack((apoly,nanbuf,aply,lastpnt))
    else:
        poly = polylist[0]
        #poly = np.vstack((poly,poly[0,:]))
        #poly = np.vstack((poly,polylist[-1][0,:]))
        apoly=skimage.measure.approximate_polygon(poly,1)
    #properties
    apoly=poly
    equivradius=blob.perimeter/(2*np.pi)
    circ_area = np.pi*equivradius**2
    Circularity = blob.area/circ_area
    cProps=(blob.area, blob.perimeter, Circularity)

    # geolocate convex hull vertices
    contour=np.zeros((apoly.shape))
    for v in range(apoly.shape[0]):
        y, x = apoly[v]
        contour[v,:] = gpx.xp2lon(x), gpx.yp2lat(y)

    contour = np.vstack((contour,contour[0,:]))
    pixcontour = np.vstack((apoly[:,1],apoly[:,0])).T

    return cProps, contour, pixcontour

def DrawContourAngles(blobs,gpx,m=plt):
    cols=['r','b','c','m','g','r','b','c','m','g']
    for rp in xrange(5):cols.extend(cols)
    for i in blobs.keys():
        cl=cols[i]
        cb=blobs[i]
        cx,cy = cb.centroid[1], cb.centroid[0]
        ### Draw angle arrow for blob
        ex = np.cos(np.deg2rad(cb.degs))*6
        ey = np.sin(np.deg2rad(cb.degs))*6
        if cb.degs<0:ey=-ey
        if hasattr(m,'drawgreatcircle'):
            cx,cy = gpx.xp2lon(cx),gpx.yp2lat(cy)
            ex,ey = gpx.xp2lon(ex),gpx.yp2lat(ey)
            mcx,mcy = m(cx,cy)
            mex,mey = m(ex,ey)
            mex2,mey2 = m(-ex,-ey)
        else:
            cx,cy = gpx.xp2lon(cx),gpx.yp2lat(cy)
            ex,ey = gpx.xp2lon(ex),gpx.yp2lat(ey)
            mcx,mcy,mex,mey = cx,cy,ex,ey
            mex2,mey2 = -ex,-ey
        plt.arrow(mcx,mcy,mex,mey,width=.07,fc=cl,ec='r')
        plt.arrow(mcx,mcy,mex2,mey2,width=.07,fc=cl,ec='r')
        txt="Cloudband\nTilt: %03.0f" %(cb.degs)
        plt.text(mcx,mcy,txt,color='c',fontsize=14.,fontweight='bold') 
        ### Draw contour
        if hasattr(m,'drawgreatcircle'):
            contour = cb.convhull
            cnx,cny = m(contour[:,0],contour[:,1])
        else:
            #cnx,cny = cb.convhullpix[:,0], cb.convhullpix[:,1]
            cnx,cny = cb.convhull[:,0], cb.convhull[:,1]
        m.plot(cnx,cny,cl,lw=3.)
        m.scatter(mcx,mcy,s=300.,c='m',marker='d')
        

# MAIN FUNCTIONS
#def GetBlobs(img,threshold,subrect=False):
def GetBlobs(img,subrect=False):
    '''Info to follow'''
    if subrect:
        x1,y1,w,h=subrect
        x1,y1,w,h=subrect
        if img.ndim>2: subimg = img[x1:x1+xw,y1:y1+yh,:]
        elif img.ndim==2: subimg = img[x1:x1+xw,y1:y1+yh]
        img=subimg
    grey=img
    #grey = skimage.img_as_ubyte(skimage.color.rgb2grey(img))
    #grey = skimage.img_as_bool( (grey>threshold[0]) & (grey<threshold[1]) )
    #grey = (grey>threshold[0]) & (grey<threshold[1])
    #plt.figure();plt.imshow(grey);plt.colorbar(shrink=.5)
    ### LABEL BLOBS
    labelled = skimage.measure.label(grey,connectivity=2,background=0)
    labelled=labelled+1
    #plt.figure();plt.imshow(labelled)
    bloblist = skimage.measure.regionprops(labelled)
    # Include labelled blobs, centroid, angle, bounding box
    blobsdict={}
    for i,blb in enumerate(bloblist): blobsdict.update({i:blb})

    return labelled, blobsdict, grey


def MetBlobs(vrb,time,hrtime,lat,lon,varstr,sub='SA',showblobs=True,\
             interact=False):
    '''mbs, blbim = MetBlobs(vrb,time,hrtime,lat,lon,varstr,sub='SA',
                             showblobs=True)
    Main blobbing loop
    USAGE: vrb (array, dim: time x lat x lon) NOTE: use plt.cm.gray_r
    therefore, whitest values are lowest, sometimes need to times vrb by (-1)'''
    plt.close('all')
    wt=100
    ### CREATE BASEMAP INSTANCES
    if showblobs: m, mfig = SAfrBasemap(lat,lon,drawstuff=True,prj='cyl')
    ###GET APROPRIATE THRESHOLDS ARE DOMAINS FROM DICTIONARY AND ASSIGN COLORMAP
    dct=filters.blobfilters[sub+'cloudband'][varstr]
    domain=dct['ROI']
    dthresh, highlow = dct['thresh']
    #data, thrs=Stretch(vrb,dthresh,tval=dstretch)
    if highlow=='low':
        data = vrb < dthresh
        mycmap=plt.cm.gray_r;#print "Reverse Gray Colormap"
    elif highlow=='high':
        data = vrb > dthresh
        mycmap=plt.cm.gray;print "Gray Colormap"
    data = skimage.img_as_ubyte(data)
    ### Initialize variables to store output
    blobslist, blobimages, chlist = [], [], []
    ### Below: allow space for 4 blobs per day, a "guestimate" initialisation
    metblobs = np.zeros((len(time)*10,14),dtype=np.float32) 
    blobtime = np.zeros((len(time)*10,4),dtype=np.int)

    if showblobs:
        bfig=plt.figure(num='Blobs')
        #bafig=plt.figure(num='BlobsAngles')
        plt.show();plt.ion()
        plt.pause(0.05)
        keyin=raw_input("Position windows as desired then press any key,\n \
        Press [Esc] at any time to quit...\
        #[Backspace] to go back one image...")
    #canvasrect = CanvasCutOut(mfig)
    #gpx = Geopix((canvasrect[2],canvasrect[3]),lat,lon)
    gpx = Geopix((len(lat),len(lon)),lat,lon)
    llgridtup=np.meshgrid(lon,lat)
    #plt.figure(num=mfig.number);plt.clf()
    t=0
    #for t in xrange(len(time)):
    ixmbs=0
    mtimstart=timer()
    while t <= (len(time)-1):
        wtimstart = timer()
        tm = time[t, :]
        hr = hrtime[t]
        humandate = "%d-%02d-%02d %02d:00" % (tm[0], tm[1], tm[2], tm[3])
        #plt.figure(num=mfig.number);plt.clf()
        #datafig=m.transform_scalar(data[t,::-1,:],lon,lat[::-1],\
        #                           len(lon),len(lat))
        #plt.figure(num=mfig.number)
        #m.imshow(datafig,mycmap,interpolation='nearest')
        #plt.clim(dstretch[0],dstretch[1])
        #plt.clim(0,255);#plt.colorbar()
        #img = fig2img(mfig,canvasrect)
        img = data[t,:,:]
        img4blob, pixss = gpx.imgsubdomainmasked(img,domain)

        #plt.figure(num=mfig.number);plt.clf()
        blbim, blobs, grey = GetBlobs(img4blob)
        if len(blobs)==0:
            print humandate,": No Blobs detected"
            if showblobs:
                plt.figure(num=mfig.number);plt.clf()
                plt.figure(num=bfig.number);plt.clf()
                if interact:
                    plt.pause(0.05)
                    d=raw_input('Press: x to stop; b to go backwards')
                else: d='nada'
                if d=='x':
                    break
                elif d=='b':
                    t=t-2
                    if t<0:
                        t=0
                t=t+1
                continue

        #FilterBlobs(dct,blobs,img,gpx,thrs)
        #if showblobs: plt.figure(num=3);plt.clf()
        FilterBlobs(dct,blobs,img,gpx)
        if len(blobs) > 0:
            #blobslist.append(blobs);blobimages.append(blbim)
            print "%s: %d CANDIDATE %s IN %s"\
                   %(humandate,len(blobs),dct['text'].upper(),varstr.upper())
            for b in blobs.keys():
                cb=blobs[b]
                miny,minx,maxy,maxx = cb.bbox
                cprops,convhull,convhullpix = BlobContour(blbim,cb,gpx)
                cb.convhull = convhull
                cb.convhullpix = convhullpix
                chlist.append(convhull)
                vrbmean,vrbmin,vrbmax=mmm_ch_val(vrb[t,:,:],llgridtup,convhull)
                metblobs[ixmbs,:]=hr,cb.label,cb.degs,\
                gpx.xp2lon(cb.centroid[1]), gpx.yp2lat(cb.centroid[0]),\
                gpx.xp2lon(minx), gpx.xp2lon(maxx), \
                gpx.yp2lat(miny), gpx.yp2lat(maxy), \
                cb.area, cprops[2], vrbmean, vrbmin, vrbmax
                ### Have removed these
                #cb.m10, cb.m01, cb.m11, cb.m20, cb.m02, cb.u11, cb.u20,
                #cb.u02, cb.n11, cb.n20, cb.n02, cb.p1 , cb.p2 
                blobtime[ixmbs,:] = tm
                ixmbs=ixmbs+1
                if ixmbs >= metblobs.shape[0]:
                    print "OOPS, ADDING MORE SPACE TO METBLOBS ARRAY"
                    addmore = np.zeros((len(time)/2,14),dtype=np.float32)
                    addtime = np.zeros((len(time)/2,4),dtype=np.int)
                    metblobs = np.append(metblobs, addmore, axis=0)
                    blobtime = np.append(blobtime, addtime, axis=0)
                #print "Candidate Details: Area= %d pixels ; Angle= %f ; Caddircularity= %f-s %f-c"\
                #%(cb.area, (cvb.cvAngle(cb)*360)/(2*np.pi), sprops[2], cprops[2])
        else:
            donada=1;print humandate

        if showblobs:
            plt.figure(num=mfig.number);plt.clf()
            plt.figure(num=bfig.number);plt.clf()
            pt1=(pixss[0],pixss[2]);pt3=(pixss[1],pixss[3])
            pt2=(pixss[0],pixss[3]);pt4=(pixss[1],pixss[2])
            pt5=pt1
            #img=skimage.color.gray2rgb(img)
            #for lc in [1,2,3,4]:
            #    exec('rr,cc=skimage.draw.line(pt%d[1],pt%d[0],pt%d[1],pt%d[0])'\
            #          %(lc,lc,lc+1,lc+1))
            #    img[rr,cc] = np.asarray([255,0,0],dtype=np.uint8)
            #    img[rr+1,cc+1] = np.asarray([255,0,0],dtype=np.uint8)
            #    img[rr-1,cc-1] = np.asarray([255,0,0],dtype=np.uint8)
            dlat=lat[1]-lat[0];dlon=lon[1]-lon[0]
            latplot = np.hstack((lat[0]-dlat/2.,lat+dlat/2.))
            lonplot = np.hstack((lon[0]-dlon/2.,lon+dlon/2.))
            plt.figure(num=mfig.number)
            plt.pcolormesh(lonplot,latplot,vrb[t,:,:],cmap=mycmap)
            plt.grid()
            DrawContourAngles(blobs,gpx,m=plt)
            plt.xlim(lonplot[0],lonplot[-1]);plt.ylim(latplot[-1],latplot[0])
            plt.draw()
            plt.figure(num=bfig.number)
            plt.pcolormesh(lonplot,latplot,blbim);plt.grid()
            plt.xlim(lonplot[0],lonplot[-1]);plt.ylim(latplot[-1],latplot[0])
            plt.draw()
            if interact:
                plt.pause(0.05)
                d=raw_input('Press: x to stop; b to go backwards')
            else:
                d='nada'
                plt.pause(.00001)
            if d=='x':
                break
            elif d=='b':
                t=t-2
                if t<0:
                    t=0
        t=t+1
        #print t," Time taken for while loop is:",timer()-wtimstart

    try:
        plt.close('all')
    except:
        dumdum='''Do nothing'''
    sm = metblobs.sum(1)
    ikeep = np.where(sm != 0)[0]
    print "Time taken for,",len(time)," timesteps in MetBlobs is:",(timer()-mtimstart)/60,"mins"

    return metblobs[ikeep,:], blobtime[ikeep,:], chlist

def MetBlobs_th(vrb,time,hrtime,lat,lon,varstr,thresh,sub='SA',showblobs=True,\
                interact=False):
    '''mbs, blbim = MetBlobs_th(vrb,time,hrtime,lat,lon,varstr,sub='SA',
                             showblobs=True)
    Main blobbing loop - but edited to be able to set threshold from wrapper
    USAGE: vrb (array, dim: time x lat x lon) NOTE: use plt.cm.gray_r
    therefore, whitest values are lowest, sometimes need to times vrb by (-1)'''
    plt.close('all')
    wt=100
    ### CREATE BASEMAP INSTANCES
    if showblobs: m, mfig = SAfrBasemap(lat,lon,drawstuff=True,prj='cyl')
    ###GET APROPRIATE THRESHOLDS ARE DOMAINS FROM DICTIONARY AND ASSIGN COLORMAP
    dct=filters.blobfilters[sub+'cloudband'][varstr]
    domain=dct['ROI']
    #dthresh, highlow = dct['thresh']
    dthresh=thresh
    highlow='low'
    #data, thrs=Stretch(vrb,dthresh,tval=dstretch)
    if highlow=='low':
        data = vrb < dthresh
        mycmap=plt.cm.gray_r;#print "Reverse Gray Colormap"
    elif highlow=='high':
        data = vrb > dthresh
        mycmap=plt.cm.gray;print "Gray Colormap"
    data = skimage.img_as_ubyte(data)
    ### Initialize variables to store output
    blobslist, blobimages, chlist = [], [], []
    ### Below: allow space for 4 blobs per day, a "guestimate" initialisation
    metblobs = np.zeros((len(time)*10,14),dtype=np.float32)
    blobtime = np.zeros((len(time)*10,4),dtype=np.int)

    if showblobs:
        bfig=plt.figure(num='Blobs')
        #bafig=plt.figure(num='BlobsAngles')
        plt.show();plt.ion()
        plt.pause(0.05)
        keyin=raw_input("Position windows as desired then press any key,\n \
        Press [Esc] at any time to quit...\
        #[Backspace] to go back one image...")
    #canvasrect = CanvasCutOut(mfig)
    #gpx = Geopix((canvasrect[2],canvasrect[3]),lat,lon)
    gpx = Geopix((len(lat),len(lon)),lat,lon)
    llgridtup=np.meshgrid(lon,lat)
    #plt.figure(num=mfig.number);plt.clf()
    t=0
    #for t in xrange(len(time)):
    ixmbs=0
    mtimstart=timer()
    while t <= (len(time)-1):
        wtimstart = timer()
        tm = time[t, :]
        hr = hrtime[t]
        humandate = "%d-%02d-%02d %02d:00" % (tm[0], tm[1], tm[2], tm[3])
        #plt.figure(num=mfig.number);plt.clf()
        #datafig=m.transform_scalar(data[t,::-1,:],lon,lat[::-1],\
        #                           len(lon),len(lat))
        #plt.figure(num=mfig.number)
        #m.imshow(datafig,mycmap,interpolation='nearest')
        #plt.clim(dstretch[0],dstretch[1])
        #plt.clim(0,255);#plt.colorbar()
        #img = fig2img(mfig,canvasrect)
        img = data[t,:,:]
        img4blob, pixss = gpx.imgsubdomainmasked(img,domain)

        #plt.figure(num=mfig.number);plt.clf()
        blbim, blobs, grey = GetBlobs(img4blob)
        if len(blobs)==0:
            print humandate,": No Blobs detected"
            if showblobs:
                plt.figure(num=mfig.number);plt.clf()
                plt.figure(num=bfig.number);plt.clf()
                if interact:
                    plt.pause(0.05)
                    d=raw_input('Press: x to stop; b to go backwards')
                else: d='nada'
                if d=='x':
                    break
                elif d=='b':
                    t=t-2
                    if t<0:
                        t=0
                t=t+1
                continue

        #FilterBlobs(dct,blobs,img,gpx,thrs)
        #if showblobs: plt.figure(num=3);plt.clf()
        FilterBlobs(dct,blobs,img,gpx)
        if len(blobs) > 0:
            #blobslist.append(blobs);blobimages.append(blbim)
            print "%s: %d CANDIDATE %s IN %s"\
                   %(humandate,len(blobs),dct['text'].upper(),varstr.upper())
            for b in blobs.keys():
                cb=blobs[b]
                miny,minx,maxy,maxx = cb.bbox
                cprops,convhull,convhullpix = BlobContour(blbim,cb,gpx)
                cb.convhull = convhull
                cb.convhullpix = convhullpix
                chlist.append(convhull)
                vrbmean,vrbmin,vrbmax=mmm_ch_val(vrb[t,:,:],llgridtup,convhull)
                metblobs[ixmbs,:]=hr,cb.label,cb.degs,\
                gpx.xp2lon(cb.centroid[1]), gpx.yp2lat(cb.centroid[0]),\
                gpx.xp2lon(minx), gpx.xp2lon(maxx), \
                gpx.yp2lat(miny), gpx.yp2lat(maxy), \
                cb.area, cprops[2], vrbmean, vrbmin, vrbmax
                ### Have removed these
                #cb.m10, cb.m01, cb.m11, cb.m20, cb.m02, cb.u11, cb.u20,
                #cb.u02, cb.n11, cb.n20, cb.n02, cb.p1 , cb.p2
                blobtime[ixmbs,:] = tm
                ixmbs=ixmbs+1
                if ixmbs >= metblobs.shape[0]:
                    print "OOPS, ADDING MORE SPACE TO METBLOBS ARRAY"
                    addmore = np.zeros((len(time)/2,14),dtype=np.float32)
                    addtime = np.zeros((len(time)/2,4),dtype=np.int)
                    metblobs = np.append(metblobs, addmore, axis=0)
                    blobtime = np.append(blobtime, addtime, axis=0)
                #print "Candidate Details: Area= %d pixels ; Angle= %f ; Caddircularity= %f-s %f-c"\
                #%(cb.area, (cvb.cvAngle(cb)*360)/(2*np.pi), sprops[2], cprops[2])
        else:
            donada=1;print humandate

        if showblobs:
            plt.figure(num=mfig.number);plt.clf()
            plt.figure(num=bfig.number);plt.clf()
            pt1=(pixss[0],pixss[2]);pt3=(pixss[1],pixss[3])
            pt2=(pixss[0],pixss[3]);pt4=(pixss[1],pixss[2])
            pt5=pt1
            #img=skimage.color.gray2rgb(img)
            #for lc in [1,2,3,4]:
            #    exec('rr,cc=skimage.draw.line(pt%d[1],pt%d[0],pt%d[1],pt%d[0])'\
            #          %(lc,lc,lc+1,lc+1))
            #    img[rr,cc] = np.asarray([255,0,0],dtype=np.uint8)
            #    img[rr+1,cc+1] = np.asarray([255,0,0],dtype=np.uint8)
            #    img[rr-1,cc-1] = np.asarray([255,0,0],dtype=np.uint8)
            dlat=lat[1]-lat[0];dlon=lon[1]-lon[0]
            latplot = np.hstack((lat[0]-dlat/2.,lat+dlat/2.))
            lonplot = np.hstack((lon[0]-dlon/2.,lon+dlon/2.))
            plt.figure(num=mfig.number)
            plt.pcolormesh(lonplot,latplot,vrb[t,:,:],cmap=mycmap)
            plt.grid()
            DrawContourAngles(blobs,gpx,m=plt)
            plt.xlim(lonplot[0],lonplot[-1]);plt.ylim(latplot[-1],latplot[0])
            plt.draw()
            plt.figure(num=bfig.number)
            plt.pcolormesh(lonplot,latplot,blbim);plt.grid()
            plt.xlim(lonplot[0],lonplot[-1]);plt.ylim(latplot[-1],latplot[0])
            plt.draw()
            if interact:
                plt.pause(0.05)
                d=raw_input('Press: x to stop; b to go backwards')
            else:
                d='nada'
                plt.pause(.00001)
            if d=='x':
                break
            elif d=='b':
                t=t-2
                if t<0:
                    t=0
        t=t+1
        #print t," Time taken for while loop is:",timer()-wtimstart

    try:
        plt.close('all')
    except:
        dumdum='''Do nothing'''
    sm = metblobs.sum(1)
    ikeep = np.where(sm != 0)[0]
    print "Time taken for,",len(time)," timesteps in MetBlobs is:",(timer()-mtimstart)/60,"mins"

    return metblobs[ikeep,:], blobtime[ikeep,:], chlist


# SPATIAL POSITIONING INFORMATION

def relpos((mbs0, mbt0), (mbs1,mbt1),twindow=13):
    '''Function calculates position vectors of second pair of metblobs relative to first pair'''
    print "Calculating relative position vectors..."
    x0, y0 = mbs0[:,3], mbs0[:,4]
    x1, y1 = mbs1[:,3], mbs1[:,4]
    tref = mbs0[:,0]
    t1 = mbs1[:,0]
    #to do: turn these vals into hours


    lenx = len(x0)
    if lenx < len(x1): lenx=len(x1)
    precalc = np.zeros((lenx,4))
    tmatch = np.zeros((lenx,1))
    ixdict={}
    for i in xrange(len(tref)):
        t=tref[i]
        tdiff = t1 - t
        ix = np.where(np.abs(tdiff) <= twindow)[0]
        nods = len(ix)
        ixdict[i]={}
        for hrdiff in np.unique(tdiff[ix]):
            ihrs = np.where(tdiff[ix]==hrdiff)[0]
            ixdict[i][hrdiff] = ix[ihrs]
        tmatch[i:i+nods] = t
        precalc[i:i+nods,0:2] = x0[i], y0[i]
        precalc[i:i+nods,2], precalc[i:i+nods,3] = x1[ix], y1[ix]

    ikeep = np.where(precalc.sum(1) != 0)[0]
    precalc=precalc[ikeep,:];tmatch=tmatch[ikeep]
    dx, dy = precalc[:,2] - precalc[:,0], precalc[:,3] - precalc[:,1]

    def showme():
        plt.figure();
        plt.hexbin(dx,dy,gridsize=(33,25),mincnt=1);plt.colorbar()
    print "...done!"
    return np.hstack((tmatch,dx[:,np.newaxis],dy[:,np.newaxis])), ixdict, showme

# TESTING AND FINE TUNING FUNCTIONS
def gethists(vrb,time,lat,lon,varstr,sub='SA',b=50,interact=False,figd=False):
    '''showhist = gethists(vrb,time,lat,lon,varstr,b=50,interact=False)

    This function produces histograms of data values at
     various stages of the MetBlobs process

     If run with interact=False, returns a callable function pointer for display
    '''
    dset, vrst, levsel, deriv = varstr.split('-')
    try:
        plt.close('all')
    except:
        dumdum='''Do nothing'''
    # Clean data up of plt.hist that doesn't handle nans
    vrbh=np.nan_to_num(vrb)

    # Get stretch intervals from dictionary
    dct=filters.blobfilters[sub+'cloudband'][varstr];domain=dct['ROI'];
    dstretch=dct['stretch'];dthresh, highlow = dct['thresh']

    #plt.figure(2);plt.hist(vrb.ravel(),bins=b)f;plt.title('Raw Data Histogram')
    #plt.plot(dthresh,0,'rv',markersize=30)
    #data, ths = Stretch(vrb,dthresh,tval=dstretch)
    #datah=np.nan_to_num(data)
    #plt.figure(3);plt.hist(data.ravel(),bins=b);plt.title('Stretched Data Histogram')
    #plt.plot(thrs,0,'rv',markersize=30)
    #if highlow=='low':
    #    mycmap=plt.cm.gray_r;print "Reverse Gray Colormap"
    #    thrs=255-ths
    #elif highlow=='high':
    #    mycmap=plt.cm.gray;print "Gray Colormap"
    #    thrs=ths

    #m, fig = SAfrBasemap(lat,lon,drawstuff=True,prj='cyl')
    #canvasrect = CanvasCutOut(fig);plt.clf()
    #
    #nsamples = 10;cnt = 0
    #w, h = canvasrect[2:4]
    #imgdata = np.zeros((nsamples,w,h))
    #ixr = randindgen(len(time),nsamples)
    #for i in ixr:
    #    #print cnt, time[i]
    #    datafig=m.transform_scalar(data[i,::-1,:],lon,lat[::-1],\
    #                               len(lon),len(lat))
    #    m.imshow(datafig,mycmap,interpolation='nearest')
    #    plt.clim(0,255)

    #    img = fig2img(fig,canvasrect)
    #    greyimg = skimage.img_as_ubyte(skimage.color.rgb2grey(img))

    #    imgdata[cnt,:,:] = np.asarray(greyimg)
    #    cnt = cnt+1;plt.clf()

    ##plt.figure(4);plt.hist(imgdata.ravel(),bins=b);plt.title('Plotted Data Histogram');plt.xlim(0,255)
    #plt.close(1)

    plt.figure(5)
    #plt.subplot(311);#plt.title(varstr.upper()+'  Hists & Threshs:')
    plt.hist(vrbh.ravel(),bins=b,normed=True,color='0.5',ec='k')
    plt.ylabel('Raw Data')
    plt.plot(dthresh,0,'k^',markersize=30)
    plt.xlim(dstretch);
    if vrst=='olr': plt.yticks(np.arange(0.002,0.016,0.004))

    #plt.subplot(312)
    #plt.hist(datah.ravel(),bins=b,normed=True,color='0.5',ec='k')
    ##plt.ylabel('Stretched Data')
    #plt.plot(ths,0,'k^',markersize=30)
    #plt.xlim((0,255));#plt.yticks(np.arange(0.002,0.016,0.004))


    #plt.subplot(313)
    #plt.hist(imgdata.ravel(),bins=b,normed=True,color='0.5',ec='k');#plt.ylabel('Plotted Data')
    #plt.plot(thrs,0,'k^',markersize=30)
    #plt.xlim((0,255));#plt.yticks(np.arange(0.002,0.016,0.004))

    showhist="what pointer?"
    if not figd: figd = os.getcwd()
    if not interact:
        figname=figd+varstr+"_valuespdf.png"
        plt.savefig(figname,dpi=150)
        #def showhist():
            #print "Try display '"+figname+"' with eye of gnome..."
            #fail=os.system("eog "+figname+" &")

    else:
        plt.show()

    return showhist
