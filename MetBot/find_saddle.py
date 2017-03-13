import numpy as np
import matplotlib.pyplot as plt
def find_saddle(data,first_guess=240,nbins=50,xlims=[100,320],method='fmin',\
                addtests=True,showplot=False,figd=False):
    '''Function to find the saddle between to peaks in a bimodal distribution.
       Usage: Best to plot a histogram of data first in order to make a good
              first guess
              The implemented deldxmin method has OLR typical proximity
              thresholds relative to first_guess value.
    '''
    import scipy.interpolate as spi
    data = np.nan_to_num(data)
    #if showplot:
    #    f,[ax1,ax2] = plt.subplots(2,1)
    #    plt.axes(ax1)
    plt.figure(num='hist')
    hy,hx,hbars=plt.hist(data.ravel(),bins=nbins,normed=True,color='0.5',ec='k')
    plt.xlim(xlims)
    plt.ylabel('Raw Data')
    plt.yticks(np.arange(0.002,0.016,0.004))
    Interpolator = spi.interp1d(hx[1:],hy,kind='cubic')
    if method=='deldxmin':
        hires_x = np.linspace(xlims[0],xlims[1],2000)
        y = Interpolator(hires_x)
        dx= y[1:]-y[:-1]
        ddx = dx[1:] - dx[:-1]
        if showplot:
            plt.axes(ax2)
            plt.plot(hires_x[:-1],dx,'b')
            plt.plot(hires_x[:-1],np.abs(dx),'g')
            plt.plot([hires_x[0],hires_x[-1]],[0,0],'b--')
            plt.twinx()
            plt.plot(hires_x[1:-1],ddx,'r')
            plt.plot([hires_x[0],hires_x[-1]],[0,0],'r--')
            plt.axes(ax1)
        iddx_mask = ddx <= 0
        data1,data2 = first_guess-30, first_guess+30
        proximity_mask= (hires_x[1:-1]<data1) | (hires_x[1:-1]>data2)
        msk = iddx_mask | proximity_mask
        idx_masked = np.ma.MaskedArray(np.abs(dx[1:]),mask=msk)
        idx_zero = idx_masked.argmin()
        olr_saddle = hires_x[idx_zero]
    elif method=='fmin':
        import scipy.optimize as spo
        [olr_saddle] = spo.fmin(Interpolator,first_guess,full_output=0,disp=0)
        olr_saddle=round(olr_saddle)
    if addtests:
        lower=olr_saddle-5
        upper=olr_saddle+5
    if showplot:
        plt.plot(olr_saddle,0,'k^',markersize=30)
        if addtests:
            plt.plot((lower,lower),(0,0.014),'k',lw=2)
            plt.plot((upper,upper),(0,0.014),'k',lw=2)
        figname=figd+'olrhist_autothresh_'+str(olr_saddle)+'.png'
        plt.savefig(figname,dpi=150)
        plt.close()
    else:
        plt.close()
    return olr_saddle

### Example usage:
if __name__=='__main__':
    olrthresh=find_saddle(olr,method='deldxmin')
    olrthresh=find_saddle(olr,method='fmin')
    olrthresh==find_saddle(olr,method='fmin',showplot=False)


