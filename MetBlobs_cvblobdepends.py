#
''' MetBlobs.py Module: Rudimentary blob detection based feature identification system
cvblob class has attributes: minx, maxx, miny, centroid(x,y), area, moments: m10, m01, m11, m20, m02, u11, u20, u02, n11, n20, n02, p1, p2
a metblobs record (mbs) has first dim time and then 27 values as below
HOURtime, Label, Angle, cX, cY, minX, maxX, minY, maxY, area, cb.m10, cb.m01, cb.m11, cb.m20, cb.m02, cb.u11, cb.u20, cb.u02, cb.n11, cb.n20, cb.n02, cb.p1, cb.p2, Circularity, varmean_in_convhull, vrbmn_in_convhull, vrbmx_in_convhull
'''
import Image
import cvblob as cvb
import os
import cv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.nxutils as nx
from matplotlib.nxutils import points_inside_poly
import mpl_toolkits.basemap as bm
from mpl_toolkits.basemap import num2date
import mytools as my
from time import time as timer
import cPickle
debugfilters=False
# SOME DEFAULT PARAMETERS
# DICTIONARY OF FILTER VALUES
from filterdict import blobfilters

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
    mask = points_inside_poly(xypts,ch)
    maskgrid = np.reshape(mask, lls)
    maskgrid = np.where(maskgrid,False,True)
    vrbma=np.ma.MaskedArray(vrb,mask=maskgrid)

    return vrbma.mean(), vrbma.min(), vrbma.max()



class Geopix:
    '''Geopix((npix_x,npix_y),lat,lon)
    Convert between image pixels and geolocations
    Methods:
    lat2pix(deg)
    lat2pix(deg)
    imgsubdomain(img,domain)
    xp2lon(xpix)
    yp2lat(ypix)'''

    def __init__(self,(npix_x,npix_y),lat,lon):
        self.rnglon = lon[-1]-lon[0]
        self.rnglat = lat[-1]-lat[0]
        self.lon0=lon[0]
        self.lat0=lat[0]
        self.lat=lat;self.lon=lon
        self.nx=npix_x;self.ny=npix_y

    def lat2pix(self,deg):
        #ideg=np.where(np.abs(self.lat-deg)==np.abs((self.lat-deg)).min())[0][0]
        #closedeg = self.lat[ideg]
        pix=(deg-self.lat0)*self.ny/self.rnglat
        #print pix
        return np.int(pix)

    def lon2pix(self,deg):
        #ideg=np.where(np.abs(self.lon-deg)==np.abs((self.lon-deg)).min())[0][0]
        #closedeg = self.lon[ideg]
        pix=(deg-self.lon0)*self.nx/self.rnglon
        #print 'lon2pix',pix,self.lon[ideg]
        return np.int(pix)

    def imgsubdomain(self,img,domain):
        if not domain:
            lon1,lon2,lat1,lat2 = self.lon[0], self.lon[-1], self.lat[0], self.lat[-1]
        else:
            lon1,lon2,lat1,lat2 = domain
        x1,x2,y1,y2 = self.lon2pix(lon1),self.lon2pix(lon2),self.lat2pix(lat1),self.lat2pix(lat2)
        subrect=(x1,y1,x2-x1,y2-y1)
        x1,y1,w,h=subrect
        subimg = cv.CreateImage((w,h), cv.IPL_DEPTH_8U, 3)
        sub = cv.GetSubRect(img,subrect)
        cv.SetData(subimg, sub.tostring())
        self.subdomain=(x1,x2,y1,y2,w,h)

        return subimg, (x1,x2,y1,y2,w,h)

    def xp2lon(self,xpix):
        x1 = self.subdomain[0]
        pix = xpix+x1
        ln = pix*self.rnglon/self.nx + self.lon0
        return np.around(ln,decimals=1)

    def yp2lat(self,ypix):
        y1 = self.subdomain[2]
        pix = ypix+y1
        lt = pix*self.rnglat/self.ny + self.lat0
        return np.around(lt,decimals=1)

def fig2img(fig,subrect):
    '''
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    Courtesy: http://www.icare.univ-lille1.fr/wiki/index.php/How_to_convert_a_matplotlib_figure_to_a_numpy_array_or_a_PIL_image'''
    # draw the renderer
    fig.canvas.draw()
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
    buf.shape = (w, h, 3); w, h, d = buf.shape
    img = Image.fromstring("RGB",(w,h),buf.tostring())
    cv_im = cv.CreateImageHeader(img.size, cv.IPL_DEPTH_8U, 3)
    cv.SetData(cv_im, img.tostring())
    x1,y1,w,h=subrect
    plotarea = cv.CreateImage((w,h), cv.IPL_DEPTH_8U, 3)
    sub = cv.GetSubRect(cv_im,subrect)
    cv.SetData(plotarea, sub.tostring())

    return plotarea

def CanvasCutOut(fig):
    fig.set_facecolor('k')
    fig.canvas.draw()
    w,h = fig.canvas.get_width_height()
    buf = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
    buf.shape = (w, h, 3)
    w, h, d = buf.shape
    img = Image.fromstring("RGB",(w,h),buf.tostring())
    cv_im = cv.CreateImageHeader(img.size, cv.IPL_DEPTH_8U, 3)
    cv.SetData(cv_im, img.tostring())
    blbim, blobs, grey = GetBlobs(cv_im,(100,255))
    gb=cvb.cvGreaterBlob(blobs);gb=blobs[gb]
    x1,y1,xw,yh=gb.minx,gb.miny,gb.maxx-gb.minx,gb.maxy-gb.miny
    subrect=(gb.minx,gb.miny,gb.maxx-gb.minx,gb.maxy-gb.miny)

    return subrect

def Stretch(arr,dthresh,fn="Linear",tval=None):
    '''stretched = Stretch(arr,fn="Linear",tval=None)
    tval = (threshmin,thresmax) stretches between user defined limits
    default tval=None stretches between min&max of array'''
    if tval:
        #print "Stretching data to 0 255..."
        arrIn=arr.copy()
        arrOut = np.where(arrIn < tval[0],tval[0],arrIn)
        arrOut = np.where(arrIn > tval[1],tval[1],arrOut)
    arrOut=arrOut-tval[0]
    m = 255.0/(tval[1]-tval[0])
    greythresh=(dthresh-tval[0])*m
    return np.round(arrOut*m), np.round(greythresh)

def Threshold(arr,tval):
    '''Black-White thresholding'''
    arrIn=arr.copy()
    arrOut = np.where(arrIn < tval,0,arrIn)
    arrOut = np.where(arrIn > tval,255,arrOut)
    return arr

def Arr2Im8UC1(arr):
    '''im = Arr2Im8UC1(arr)

    Converts ndarray slice into 8-bit grayscale image.
    Needed to apply cv.ExtractSURF to image

    Input: arr is a 2-d ndarray
    Returns: im an CvMat(type=CV_8UC1)'''
    ex = arr.copy()
    mat = cv.fromarray(ex)
    im = cv.CreateMat(mat.rows,mat.cols,cv.CV_8UC1)
    cv.Convert(mat,im)

    return im

def Load2IPL8U(imfile,thrs):
    img = cv.LoadImage(imfile, 1)
    #convert to a greyscale image, and set intensity thresholds
    grey = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_8U, 1)
    cv.CvtColor(img, grey, cv.CV_BGR2GRAY)
    cv.Threshold(grey, grey, thrs[0], thrs[1], cv.CV_THRESH_BINARY)

    return img, grey

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

def FilterByAngle(dct,blobs,img,gpx,thrs):
    minAngle, maxAngle = dct['angle']
    if 'angleROI' in dct.keys():   # Use an angle from blobs in subdomain of one used to get blobs passed to FilterByAngle
        angledomain = dct['angleROI']
        bigsubdomain = gpx.subdomain
        img4blob, pixss = gpx.imgsubdomain(img,angledomain)
        smallsubdomain = gpx.subdomain
        angblbim, angblobs, anggrey = GetBlobs(img4blob,(thrs,255))
        dstold = 100000
        for i in blobs.keys():
            gpx.subdomain=bigsubdomain
            xb,yb = gpx.xp2lon(blobs[i].centroid[0]), gpx.yp2lat(blobs[i].centroid[1])
            ij = 1
            for j in angblobs.keys():    # Associate angleblobs with blobs
                gpx.subdomain=smallsubdomain
                xab,yab = gpx.xp2lon(angblobs[j].centroid[0]), gpx.yp2lat(angblobs[j].centroid[1])
                dst = np.sqrt((xab-xb)**2 + (yab-yb)**2) # Know this should be distance on sphere but just need crude approximation
                if dst < dstold: # Test to get angleblob centroid closest to blob centroid
                    ijmatch = j
            degs=(cvb.cvAngle(blobs[i])*360)/(2*np.pi)
            angdegs=(cvb.cvAngle(angblobs[ijmatch])*360)/(2*np.pi)
            if ( degs < minAngle or degs > maxAngle) & ( angdegs < minAngle or angdegs > maxAngle): # If associated angleblob angle fails, delete parent blob from dictionary
                del blobs[i]
                #print "Failed angles"
        gpx.subdomain=bigsubdomain
    else:          # Nice and simple unlike above
        for i in blobs.keys():
            degs=(cvb.cvAngle(blobs[i])*360)/(2*np.pi)
            if ( degs < minAngle or degs > maxAngle):
                del blobs[i]

def FilterByLatextent(dct,blobs,gpx):
    latN,latS=dct['latextent']
    pixN, pixS = gpx.lat2pix(latN), gpx.lat2pix(latS)
    if hasattr(gpx,'subdomain'):
        y1, y2 = gpx.subdomain[2], gpx.subdomain[3]
        if pixN < y1:
            pixN = y1;
        if pixS > y2:
            pixS = y2;

    mn, mx = pixN+5-y1, pixS-5-y1
    #print "Lat pixs:", pixN, pixS
    for i in blobs.keys():
        #print 'Current iteration',i, blobs[i].label
        if (blobs[i].miny > mn or blobs[i].maxy < mx):
            #print 'Deleting blob no.',i,' labelled:',blobs[i].label
            del blobs[i]
            #print "Failed Latextent"

def testshow(im):
    cv.StartWindowThread()
    cv.NamedWindow('test',1)
    cv.ShowImage('test',im)
    cv.WaitKey(1)

# MAPPING AND PLOTTING HELPERS
def SAfrBasemap(lat,lon,drawstuff=False,prj='cyl',fno=1,rsltn='c'):
    '''m, f = SAfrBasemap(lat,lon,drawstuff=False,prj='cyl',fno=1)

    This creates a basemap instance from lat's and lon's provided.
    Specific for this application of metblobs, so currently using proj='cyl', however
    albers equal area is here uncommented.
    USAGE: lat, lon
    RETURNS: m, basemap object pointer (handle in matlab language)
             f, pointer figure'''
    xy=(lat.min(),lat.max(),lon.min(),lon.max())
    nx = len(lon); ny = len(lat)

    if prj=='cyl':
        m = bm.Basemap(llcrnrlon=xy[2],llcrnrlat=xy[0],urcrnrlon=xy[3],urcrnrlat=xy[1],resolution=rsltn,area_thresh=10000.,projection='cyl')
    if prj=='aea':
        m = bm.Basemap(llcrnrlon=xy[2]-5.,llcrnrlat=xy[0],urcrnrlon=xy[3],urcrnrlat=xy[1]+5,resolution=rsltn,projection='aea', lat_1=-45.,lat_2=0.,lon_0=40.)

    ### SET UP FIGURE AND MERIDIANS
    delon = 10.
    meridians = np.arange(10.,360.,delon)
    delat = 5.
    circles = np.arange(0.,90.+delat,delat).tolist()+np.arange(-delat,-90.-delat,-delat).tolist()
    f1 = plt.figure(fno,figsize=[12.0, 8.0])
    #ax = f1.add_axes([0.1,0.1,0.7,0.7])
    #pos = ax.get_position()
    #l, b, w, h = pos.bounds
    #cax = plt.axes([l+w+0.075, b, 0.05, h]) # setup colorbar axes.
    #plt.colorbar(cax=cax) # draw colorbar
    #plt.axes(ax)

    if drawstuff:
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(circles,linewidth='0.1',labels=[1,0,0,0])
        m.drawmeridians(meridians,linewidth='0.1',labels=[0,0,0,1])

    return m, f1

myfont=cv.InitFont(cv.CV_FONT_HERSHEY_COMPLEX_SMALL, 1., 1., shear=0, thickness=2, lineType=8)

def BlobContour(imgOut,blob,gpx,render=True,chtype='simple'):
    '''chtype=convexhull or simple (which is more complex)'''
    polygon = cvb.cvConvertChainCodesToPolygon(blob.contour)
    sp = cvb.cvSimplifyPolygon(polygon, 10.)
    cp = cvb.cvPolygonContourConvexHull(sp)

    if render:
        cvb.cvRenderContourPolygon(sp, imgOut,cv.CV_RGB(0, 0, 255))
        cvb.cvRenderContourPolygon(cp, imgOut,cv.CV_RGB(0, 255, 0))

    #properties
    spProps=(cvb.cvContourPolygonArea(sp), cvb.cvContourPolygonPerimeter(sp), cvb.cvContourPolygonCircularity(sp))
    cpProps=(cvb.cvContourPolygonArea(cp), cvb.cvContourPolygonPerimeter(cp), cvb.cvContourPolygonCircularity(cp))

    # geolocate convex hull vertices
    if chtype=='convexhull':
        polyg=cp
    elif chtype=='simple':
        polyg=sp

    contour=np.zeros((polyg.size(),2))
    for v in range(polyg.size()):
        x, y = polyg.pop()
        contour[v,:] = gpx.xp2lon(x), gpx.yp2lat(y)

    contour = np.vstack((contour,contour[0,:]))

    return spProps, cpProps, contour


# MAIN FUNCTIONS
def GetBlobs(img,threshold,subrect=False):
    '''Info to follow'''
    if subrect:
        x1,y1,w,h=subrect
        subimg = cv.CreateImage((w,h), cv.IPL_DEPTH_8U, 1)
        sub = cv.GetSubRect(img,subrect)
        cv.SetData(subimg, sub.tostring())
        img=subimg

    grey = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_8U, 1)
    cv.CvtColor(img, grey, cv.CV_RGB2GRAY)
    cv.Threshold(grey, grey, threshold[0], threshold[1], cv.CV_THRESH_BINARY)
    ### LABEL BLOBS
    IPL_DEPTH_LABEL=32
    labelImg = cv.CreateImage(cv.GetSize(img), IPL_DEPTH_LABEL, 1)
    blobs = cvb.CvBlobs()
    result = cvb.cvLabel(grey, labelImg, blobs)
    ### RENDER LABELLED BLOBS
    imgOut = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_8U, 3)
    cv.Zero(imgOut)
    cvb.cvRenderBlobs(labelImg, blobs, img, imgOut, cvb.CV_BLOB_RENDER_COLOR|cvb.CV_BLOB_RENDER_CENTROID|cvb.CV_BLOB_RENDER_BOUNDING_BOX|cvb.CV_BLOB_RENDER_ANGLE, 1.0)
    #cv.SaveImage('tmp/'+imf, imgOut)

    return imgOut, blobs, grey


def MetBlobs(vrb,time,hrtime,lat,lon,varstr,sub='SA',showblobs=True,contype='simple'):
    '''mbs, blbim = MetBlobs(vrb,time,hrtime,lat,lon,varstr,sub='SA',showblobs=True,contype='simple')

    Main blobbing loop
    USAGE: vrb (array, dim: time x lat x lon) NOTE: use plt.cm.gray_r therefore,
    whitest values are lowest, sometimes need to times vrb by (-1)'''
    plt.close('all')
    wt=100
    ### CREATE BASEMAP INSTANCES
    m, fig = SAfrBasemap(lat,lon,drawstuff=True,prj='cyl')
    canvasrect = CanvasCutOut(fig);
    if False:
        plt.show();return

    gpx = Geopix((canvasrect[2],canvasrect[3]),lat,lon)
    llgridtup=np.meshgrid(lon,lat)
    ### GET APROPRIATE THRESHOLDS ARE DOMAINS FROM DICTIONARY AND ASSIGN COLORMAP
    dct=blobfilters[sub+'cloudband'][varstr]
    domain=dct['ROI']
    dstretch=dct['stretch']
    dthresh, highlow = dct['thresh']
    data, thrs=Stretch(vrb,dthresh,tval=dstretch)
    if highlow=='low':
        mycmap=plt.cm.gray_r;#print "Reverse Gray Colormap"
        thrs=255-thrs
    elif highlow=='high':
        mycmap=plt.cm.gray;print "Gray Colormap"

    blobslist=[]
    blobimages=[]
    chlist=[]
    metblobs = np.zeros((len(time)*10,27),dtype=np.float32) # this essentially allows space for 4 blobs per day, a "guestimate" initialisation
    blobtime = np.zeros((len(time)*10,4),dtype=np.int)

    if showblobs:
        cv.StartWindowThread()
        cv.NamedWindow("MetBlobs",1);cv.MoveWindow("MetBlobs",0,canvasrect[3]+30)
        cv.NamedWindow("Raw OLR",1);cv.MoveWindow("Raw OLR",0,0)
        print "Position windows as desired then press any key,\n \
        Press [Esc] at any time to quit...\
        #[Backspace] to go back one image...";cv.WaitKey(1)
    plt.clf()
    t=0
    #for t in xrange(len(time)):
    ixmbs=0
    mtimstart=timer()
    while t <= (len(time)-1):
        wtimstart=timer()
        tm=time[t,:]
        hr=hrtime[t]
        #humandate=num2date(time[t],units="hours since 1-1-1 00:00:0.0",calendar='gregorian')
        humandate="%d-%02d-%02d %02d:00" %(tm[0],tm[1],tm[2],tm[3])
        #print num2date(time[t],units="hours since 1-1-1 00:00:0.0",calendar='gregorian')
        datafig=m.transform_scalar(data[t,::-1,:],lon,lat[::-1],len(lon),len(lat))
        m.imshow(datafig,mycmap,interpolation='nearest')
        #plt.clim(dstretch[0],dstretch[1])
        plt.clim(0,255)
        img = fig2img(fig,canvasrect)
        img4blob, pixss = gpx.imgsubdomain(img,domain)

        ### OR ###
        #x, y = m(*np.meshgrid(lon,lat))
        #mp=m.pcolormesh(x,y,olr[t,:,:],cmap=mycmap)
        #img = blb.fig2img(f1)
        #pim=img.crop((139,184,700,540))

        #img.show()         # this is useful when needing to debug

        fig.clf()

        blbim, blobs, grey = GetBlobs(img4blob,(thrs,255))
        if len(blobs)==0:
            print humandate,": No Blobs detected"
            if showblobs:
                cv.ShowImage("Raw OLR",img)
                d=cv.WaitKey(wt)
                if d==27 or d==1048603:  # ESC has second value when numlock is on for my keyboard
                    break
                elif d==1113864 or d==10:
                    t=t-3
                    if t<0:
                        t=0
                t=t+1
                continue

        if False:
            if len(blobs) > 0:
                gb=blobs[cvb.cvGreaterBlob(blobs)];sprops, cprops,convhull = BlobContour(blbim,gb,gpx,chtype=contype)
                #print "Greater Blob: Area= %d pixels ; Angle= %f ; Circularity= %f-simple %f-complex"\
                #%(gb.area, (cvb.cvAngle(gb)*360)/(2*np.pi), sprops[2], cprops[2])

        FilterBlobs(dct,blobs,img,gpx,thrs)
        if len(blobs) > 0:
            #blobslist.append(blobs);blobimages.append(blbim)
            print humandate,": ",len(blobs)," CANDIDATE ",dct['text'].upper()," IN ",varstr.upper()
            for b in blobs.keys():
                cb=blobs[b]
                sprops, cprops, convhull = BlobContour(blbim,cb,gpx,chtype=contype)
                chlist.append(convhull)
                vrbmean, vrbmin, vrbmax = mmm_ch_val(vrb[t,:,:],llgridtup,convhull)
                metblobs[ixmbs,:] = hr, cb.label, (cvb.cvAngle(cb)*360.0)/(2*np.pi), gpx.xp2lon(cb.centroid[0]), gpx.yp2lat(cb.centroid[1]), gpx.xp2lon(cb.minx), gpx.xp2lon(cb.maxx), gpx.yp2lat(cb.miny), gpx.yp2lat(cb.maxy), cb.area, cb.m10, cb.m01, cb.m11, cb.m20, cb.m02, cb.u11, cb.u20, cb.u02, cb.n11, cb.n20, cb.n02, cb.p1, cb.p2, cprops[2], vrbmean, vrbmin, vrbmax
                blobtime[ixmbs,:] = tm
                ixmbs=ixmbs+1
                if ixmbs >= metblobs.shape[0]:
                    print "OOPS, ADDING MORE SPACE TO METBLOBS ARRAY"
                    addmore = np.zeros((len(time)/2,27),dtype=np.float32)
                    addtime = np.zeros((len(time)/2,4),dtype=np.int)
                    metblobs = np.append(metblobs, addmore, axis=0)
                    blobtime = np.append(blobtime, addtime, axis=0)
                #print "Candidate Details: Area= %d pixels ; Angle= %f ; Caddircularity= %f-s %f-c"\
                #%(cb.area, (cvb.cvAngle(cb)*360)/(2*np.pi), sprops[2], cprops[2])
                cv.PutText(img,dct['text'],(int(cb.centroid[0])-50+pixss[0],int(cb.centroid[1])+pixss[2]),myfont,4)
        else:
            donada=1;print humandate

        if showblobs:
            pt1=(pixss[0],pixss[2]);pt2=(pixss[1],pixss[3])
            cv.Rectangle(img,pt1,pt2,cv.RGB(225,0,0),thickness=2)
            cv.ShowImage("Raw OLR",img)
            cv.ShowImage("MetBlobs",blbim)
            d=cv.WaitKey(wt)
            if d==27 or d==1048603:  # ESC has second value when numlock is on for my keyboard
                break
            elif d==1113864 or d==10 or d==65288:
                t=t-2
                if t<0:
                    t=0
        t=t+1
        #print t," Time taken for while loop is:",timer()-wtimstart

    if showblobs:
        cv.DestroyAllWindows()
    try:
        plt.close('all')
    except:
        dumdum='''Do nothing'''
    sm = metblobs.sum(1)
    ikeep = np.where(sm != 0)[0]
    print "Time taken for,",len(time)," timesteps in MetBlobs is:",(timer()-mtimstart)/60,"mins"

    return metblobs[ikeep,:], blobtime[ikeep,:], chlist


def FilterBlobs(dct,blobs,img,gpx,thrs):
    '''Filters blobs by predefined criteria:
    These are contained in blobfilters (dict type)'''
    # Area Filter
    if dct['area']:
        mn,mx=dct['area']
        cvb.cvFilterByArea(blobs,mn,mx)
    # Latitutdinal Extent
    if dct['latextent']:
        FilterByLatextent(dct,blobs,gpx)
    # Angle Filter
    if dct['angle']:
        FilterByAngle(dct,blobs,img,gpx,thrs)


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
def gethists(vrb,time,lat,lon,varstr,sub='SA',b=50,interact=False):
    '''showhist = gethists(vrb,time,lat,lon,varstr,b=50,interact=False)

    This function produces histograms of data values at
     various stages of the MetBlobs process

     If run with: interact=False, returns a callable function pointer for display'''
    dset, vrst, levsel, deriv = varstr.split('-')
    try:
        plt.close('all')
    except:
        dumdum='''Do nothing'''
    # Clean data up of plt.hist that doesn't handle nans
    vrbh=np.nan_to_num(vrb)

    # Get stretch intervals from dictionary
    dct=blobfilters[sub+'cloudband'][varstr];domain=dct['ROI'];
    dstretch=dct['stretch'];dthresh, highlow = dct['thresh']

    #plt.figure(2);plt.hist(vrb.ravel(),bins=b);plt.title('Raw Data Histogram')
    #plt.plot(dthresh,0,'rv',markersize=30)
    data, ths = Stretch(vrb,dthresh,tval=dstretch)
    datah=np.nan_to_num(data)
    #plt.figure(3);plt.hist(data.ravel(),bins=b);plt.title('Stretched Data Histogram')
    #plt.plot(thrs,0,'rv',markersize=30)
    if highlow=='low':
        mycmap=plt.cm.gray_r;print "Reverse Gray Colormap"
        thrs=255-ths
    elif highlow=='high':
        mycmap=plt.cm.gray;print "Gray Colormap"
        thrs=ths

    m, fig = SAfrBasemap(lat,lon,drawstuff=True,prj='cyl');canvasrect = CanvasCutOut(fig);plt.clf()

    nsamples = 10;cnt = 0
    w, h = canvasrect[2:4]
    imgdata = np.zeros((nsamples,h,w))
    ixr = my.randindgen(len(time),nsamples)
    for i in ixr:
        #print cnt, time[i]
        datafig = m.transform_scalar(data[i,::-1,:],lon,lat[::-1],len(lon),len(lat))
        m.imshow(datafig,mycmap,interpolation='nearest')
        plt.clim(0,255)

        img = fig2img(fig,canvasrect)
        grey = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_8U, 1)
        cv.CvtColor(img, grey, cv.CV_RGB2GRAY)

        greypil = Image.fromstring("L", cv.GetSize(grey), grey.tostring())
        imgdata[cnt,:,:] = np.asarray(greypil)
        cnt = cnt+1;plt.clf()

    #plt.figure(4);plt.hist(imgdata.ravel(),bins=b);plt.title('Plotted Data Histogram');plt.xlim(0,255)
    plt.close(1)

    plt.figure(5)
    plt.subplot(311);#plt.title(varstr.upper()+'  Hists & Threshs:')
    plt.hist(vrbh.ravel(),bins=b,normed=True,color='0.5',ec='k');#plt.ylabel('Raw Data');
    plt.plot(dthresh,0,'k^',markersize=30)
    plt.xlim(dstretch);
    if vrst=='olr': plt.yticks(np.arange(0.002,0.016,0.004))

    plt.subplot(312)
    plt.hist(datah.ravel(),bins=b,normed=True,color='0.5',ec='k');#plt.ylabel('Stretched Data')
    plt.plot(ths,0,'k^',markersize=30)
    plt.xlim((0,255));#plt.yticks(np.arange(0.002,0.016,0.004))


    plt.subplot(313)
    plt.hist(imgdata.ravel(),bins=b,normed=True,color='0.5',ec='k');#plt.ylabel('Plotted Data')
    plt.plot(thrs,0,'k^',markersize=30)
    plt.xlim((0,255));#plt.yticks(np.arange(0.002,0.016,0.004))

    showhist="what pointer?"
    if not interact:
        figd='/home/neil/work/computervision/metblobs/'+sub+'/rawhists/'
        if dset=='hadam3p':figd='/home/neil/work/computervision/metblobs/'+sub+'/rawhists/flavours/'
        figname=figd+varstr+"_valuespdf.png"
        plt.savefig(figname,dpi=150)
        #def showhist():
            #print "Try display '"+figname+"' with eye of gnome..."
            #fail=os.system("eog "+figname+" &")

    else:
        plt.show()

    return showhist
