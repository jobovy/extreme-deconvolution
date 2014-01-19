import numpy
from extreme_deconvolution import extreme_deconvolution
from galpy.util import bovy_plot
def test1_ngerrors():
    #Generate data
    ndata= 10001
    ngauss= 1
    ydata= numpy.random.normal(scale=[1.,2.],size=(ndata,2))*numpy.sqrt(2.)
    ycovar= numpy.ones((ndata,2))*0.
    ngamp= numpy.ones((ndata,ngauss))/ngauss
    ngmean= numpy.zeros((ndata,ngauss,2))
    ngcovar= numpy.ones((ndata,ngauss,2))
    xamp= numpy.ones(1)/1.
    xmean= numpy.array([[0.,0.]])
    xcovar= numpy.array([[[ 0.03821028, 0.02108113],
                          [ 0.02108113,  0.03173839]]])
#    """
    xamp= numpy.ones(2)/2.
    xmean= numpy.array([[0.,0.],[1.,-1.]])
    xcovar= numpy.array([[[ 0.03821028, 0.02108113],
                          [ 0.02014796,  0.03173839]],
                         [[ 0.06219194,  0.02302473],
                          [ 0.02738021,  0.06778009]]])
#    """
    l= extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,
                             ng=True,
                             ngamp=ngamp,
                             ngmean=ngmean,
                             ngcovar=ngcovar)
    print l
    print xamp, xmean, xcovar
    
    
def test_ngerrors():
    samples= False
    if samples:
        ngauss= 10
    else:
        ngauss= 2
    #Generate data
    ndata= 10001
    ydata= numpy.random.normal(size=(ndata,1))
    #Add noise
    for ii in range(ndata):
        if not samples:
            if numpy.random.uniform() < 0.5:
                ydata[ii,0]+= numpy.random.normal()+5.
            else:
                ydata[ii,0]+= numpy.random.normal()-3.
    #bovy_plot.bovy_print()
    #bovy_plot.bovy_hist(ydata,bins=101,histtype='step',color='k')
    #bovy_plot.bovy_end_print('/Users/bovy/Desktop/test.png')
    ycovar= numpy.ones((ndata,2))*0.
    ngamp= numpy.ones((ndata,ngauss,1))/ngauss
    ngmean= numpy.zeros((ndata,ngauss,1))
    if samples:
        for ii in range(ndata):
            for jj in range(ngauss):
                if numpy.random.uniform() < 0.5:
                    ngmean[ii,jj,0]= ydata[ii,0]+(numpy.random.normal()+5.)
                else:
                    ngmean[ii,jj,0]= ydata[ii,0]+(numpy.random.normal()-3.)
            ydata[ii,0]= 0.
        ngcovar= numpy.zeros((ndata,ngauss,1))
    else:
        ngmean[:,0,0]= 5.
        ngmean[:,1,0]= -3.
        ngcovar= numpy.ones((ndata,ngauss,1))
    xamp= numpy.ones(1)/1.
    xmean= numpy.array([[0.]])
    xcovar= numpy.array([[[0.03821028]]])
    """
    xamp= numpy.ones(2)/2.
    xmean= numpy.array([[0.,0.],[1.,-1.]])
    xcovar= numpy.array([[[ 0.03821028, 0.02108113],
                          [ 0.02014796,  0.03173839]],
                         [[ 0.06219194,  0.02302473],
                          [ 0.02738021,  0.06778009]]])
    """
    l= extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,
                             ng=True,
                             ngamp=ngamp,
                             ngmean=ngmean,
                             ngcovar=ngcovar)
    print l
    print xamp, xmean, xcovar
    
    

if __name__ == '__main__':
    numpy.random.seed(1)
    test_ngerrors()
