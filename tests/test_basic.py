# test_basic.py: some basic tests of the code
import numpy
numpy.random.seed(2)
from extreme_deconvolution import extreme_deconvolution

def test_single_gauss_1d_nounc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 1001
    ydata= numpy.atleast_2d(numpy.random.normal(size=ndata)).T
    ycovar= numpy.zeros_like(ydata)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata)+numpy.std(ydata))
    initcovar= numpy.atleast_3d(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 3./numpy.sqrt(ndata)
    assert numpy.fabs(initmean-0.) < tol, 'XD does not recover correct mean for single Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar-1.) < tol, 'XD does not recover correct variance for single Gaussian w/o uncertainties'
    return None
