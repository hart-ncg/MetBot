#
import Image
import cvblob as cvb
import cv
import os
import cv
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import mytools as my

# DICTIONARY OF FILTER VALUES
blobfilters={}
blobfilters['cloudband']={'olr':\
 {'minpixels': 5000, 'maxpixels': 1000000,\
 'mindeg':15.0, 'maxdeg': 90.0,\
 'ROI': (7.5, 70,-20,-35),\
 'stretch': (80,300),\
 'thresh': 130 }}

# Helper functions
class Geopix:
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
        lon1,lon2,lat1,lat2 = domain
        x1,x2,y1,y2 = self.lon2pix(lon1),self.lon2pix(lon2),self.lat2pix(lat1),self.lat2pix(lat2)
        subrect=(x1,y1,x2-x1,y2-y1)
        x1,y1,w,h=subrect
        subimg = cv.CreateImage((w,h), cv.IPL_DEPTH_8U, 3)
        sub = cv.GetSubRect(img,subrect)
        cv.SetData(subimg, sub.tostring())
        return subimg, (x1,x2,y1,y2,w,h)

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
    gb=cvb.GreaterBlob(blobs);gb=blobs[gb]
    x1,y1,xw,yh=gb.minx,gb.miny,gb.maxx-gb.minx,gb.maxy-gb.miny
    subrect=(gb.minx,gb.miny,gb.maxx-gb.minx,gb.maxy-gb.miny)

    return subrect

def Stretch(arr,dthresh,fn="Linear",tval=None):
    '''stretched = Stretch(arr,fn="Linear",tval=None)
    tval = (threshmin,thresmax) stretches between user defined limits
    default tval=None stretches between min&max of array'''
    if tval:
        print "Stretching data to 0 255..."
        arrIn=arr.copy()
        arrOut = np.where(arrIn < tval[0],tval[0],arrIn)
        arrOut = np.where(arrIn > tval[1],tval[1],arrOut)
    arrOut=arrOut-tval[0]
    m = 255.0/(tval[1]-tval[0])
    greythresh=255*dthresh/(tval[1]-tval[0])
    return np.round(arrOut*m), greythresh

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

def FilterByAngle(blobs, minAngle, maxAngle):
    for i in blobs.keys():
        degs=(cvb.Angle(blobs[i])*360)/(2*np.pi)
        if ( degs < minAngle or degs > maxAngle):
            del blobs[i]

def testshow(im):
    cv.StartWindowThread()
    cv.NamedWindow('test',1)
    cv.ShowImage('test',im)
# MAPPING AND PLOTTING HELPERS
def SAfrBasemap(lat,lon,drawstuff=False,prj='cyl'):
    '''m, f = SAfrBasemap(lat,lon,drawstuff=False,prj='cyl')

    This creates a basemap instance from lat's and lon's provided.
    Specific for this application of metblobs, so currently using proj='cyl', however
    albers equal area is here uncommented.
    USAGE: lat, lon
    RETURNS: m, basemap object pointer (handle in matlab language)
             f, pointer figure'''
    xy=(lat.min(),lat.max(),lon.min(),lon.max())
    nx = len(lon); ny = len(lat)

    if prj=='cyl':
        m = bm.Basemap(llcrnrlon=xy[2],llcrnrlat=xy[0],urcrnrlon=xy[3],urcrnrlat=xy[1],resolution='c',area_thresh=10000.,projection='cyl')
    if prj=='aea':
        m = bm.Basemap(llcrnrlon=xy[2]-5.,llcrnrlat=xy[0],urcrnrlon=xy[3],urcrnrlat=xy[1]+5,resolution='c',projection='aea', lat_1=-45.,lat_2=-0.,lon_0=40.)

    ### SET UP FIGURE AND MERIDIANS
    delon = 10.
    meridians = np.arange(10.,360.,delon)
    delat = 5.
    circles = np.arange(0.,90.+delat,delat).tolist()+np.arange(-delat,-90.-delat,-delat).tolist()
    f1 = plt.figure(1,figsize=[12.0, 8.0])
    #ax = f1.add_axes([0.1,0.1,0.7,0.7])
    #pos = ax.get_position()
    #l, b, w, h = pos.bounds
    #cax = plt.axes([l+w+0.075, b, 0.05, h]) # setup colorbar axes.
    #plt.colorbar(cax=cax) # draw colorbar
    #plt.axes(ax)
    if drawstuff:
        m.drawcoastlines()
        m.drawcountries()
        m.drawparallels(circles,linewidth='0.1',labels=[1,1,0,0])
        m.drawmeridians(meridians,linewidth='0.1',labels=[0,0,1,1])

    return m, f1

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
    labelImg = cv.CreateImage(cv.GetSize(img), cvb.IPL_DEPTH_LABEL, 1)
    blobs = cvb.Blobs()
    result = cvb.Label(grey, labelImg, blobs)
    ### RENDER LABELLED BLOBS
    imgOut = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_8U, 3)
    cv.Zero(imgOut)
    cvb.RenderBlobs(labelImg, blobs, img, imgOut, cvb.CV_BLOB_RENDER_COLOR|cvb.CV_BLOB_RENDER_CENTROID|cvb.CV_BLOB_RENDER_BOUNDING_BOX|cvb.CV_BLOB_RENDER_ANGLE, 1.0)
    #cv.SaveImage('tmp/'+imf, imgOut)

    return imgOut, blobs, grey

def MetBlobs(vrb,time,lat,lon,varstr,showblobs=True):
    '''Main blobbing loop
    USAGE: thrs tuple'''
    ### CREATE BASEMAP INSTANCES
    m, fig = SAfrBasemap(lat,lon,drawstuff=True,prj='cyl')
    canvasrect = CanvasCutOut(fig);
    if False:
        plt.show();return

    gpx = Geopix((canvasrect[2],canvasrect[3]),lat,lon)
    ### GET APROPRIATE THRESHOLDS ARE DOMAINS FROM DICTIONARY
    dct=blobfilters['cloudband'][varstr]
    domain=dct['ROI']
    dstretch=dct['stretch']
    dthresh=dct['thresh']
    data, thrs=Stretch(vrb,dthresh,tval=dstretch)
    metblobs=[]
    blobimages=[]

    if showblobs:
        cv.StartWindowThread()
        cv.NamedWindow("MetBlobs",1)
        cv.NamedWindow("Raw OLR",1)
        cv.NamedWindow("Candidate Cloudband",1)
        print "Position windows as desired then press any key,\n \
        Press [Esc] at any time to quit...";cv.WaitKey()
    plt.clf()
    for t in xrange(len(time)):
        print time[t]
        datafig=m.transform_scalar(data[t,::-1,:],lon,lat[::-1],len(lon),len(lat))
        m.imshow(datafig,plt.cm.gray_r,interpolation='bicubic')
        plt.clim(dstretch[0],dstretch[1])
        img = fig2img(fig,canvasrect)
        img4blob, pixss = gpx.imgsubdomain(img,domain)
        ### OR ###
        #x, y = m(*np.meshgrid(lon,lat))
        #mp=m.pcolormesh(x,y,olr[t,:,:],cmap=plt.cm.gray_r)
        #img = blb.fig2img(f1)
        #pim=img.crop((139,184,700,540))

        #img.show()         # this is useful when needing to debug

        fig.clf()

        blbim, blobs, grey = GetBlobs(img4blob,(thrs,255))
        if showblobs:
            #cv.ResetImageROI
            pt1=(pixss[0],pixss[2]);pt2=(pixss[1],pixss[3])
            cv.Rectangle(img,pt1,pt2,cv.RGB(225,0,0),thickness=2)
            cv.ShowImage("Raw OLR",img)
            cv.ShowImage("MetBlobs",blbim)
            d=cv.WaitKey()
            if d==27 or d==1048603:  # ESC has second value when numlock is on for my keyboard
                break

        if len(blobs)==0:
            print "No Blobs detected"
            continue

        gb=blobs[cvb.GreaterBlob(blobs)]
        print "Greatest Blob covers: %d pixels @ angle of %f" %(gb.area, (cvb.Angle(gb)*360)/(2*np.pi))
        FilterBlobs(blobs,varstr)
        if len(blobs) > 0:
            metblobs.append(blobs)
            print "CANDIDATE CLOUDBAND IN OLR"
            if showblobs:
                cv.ShowImage("Candidate Cloudband",blbim)
                #cv.WaitKey()
    cv.DestroyAllWindows()

    return metblobs, blobimages

def FilterBlobs(blobs,varstr,wtype='cloudband'):
    '''Filters blobs by predefined criteria:
    These are contained in blobfilters (dict type)'''
    typefilters=blobfilters[wtype]
    vals=typefilters[varstr]
    mn=vals['minpixels'];mx=vals['maxpixels']
    cvb.FilterByArea(blobs,mn,mx)
    mn=vals['mindeg'];mx=vals['maxdeg']
    FilterByAngle(blobs,mn,mx)
