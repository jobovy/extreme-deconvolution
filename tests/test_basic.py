# test_basic.py: some basic tests of the code
import numpy
numpy.random.seed(2)
from extreme_deconvolution import extreme_deconvolution

def test_single_gauss_1d_nounc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
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

def test_single_gauss_1d_constunc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    ydata= numpy.atleast_2d(numpy.random.normal(size=ndata)).T
    ycovar= numpy.ones_like(ydata)*0.25
    ydata+= numpy.atleast_2d(numpy.random.normal(size=ndata)).T\
        *numpy.sqrt(ycovar)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata)+numpy.std(ydata))
    initcovar= numpy.atleast_3d(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 3./numpy.sqrt(ndata)
    assert numpy.fabs(initmean-0.) < tol, 'XD does not recover correct mean for single Gaussian w/ constant uncertainties'
    assert numpy.fabs(initcovar-1.) < tol, 'XD does not recover correct variance for single Gaussian w/ constant uncertainties'
    return None

def test_single_gauss_1d_varunc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    ydata= numpy.atleast_2d(numpy.random.normal(size=ndata)).T
    ycovar= numpy.ones_like(ydata)*\
        numpy.atleast_2d(numpy.random.uniform(size=ndata)).T
    ydata+= numpy.atleast_2d(numpy.random.normal(size=ndata)).T\
        *numpy.sqrt(ycovar)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata)+numpy.std(ydata))
    initcovar= numpy.atleast_3d(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 5./numpy.sqrt(ndata)
    assert numpy.fabs(initmean-0.) < tol, 'XD does not recover correct mean for single Gaussian w/ uncertainties'
    assert numpy.fabs(initcovar-1.) < tol, 'XD does not recover correct variance for single Gaussian w/ uncertainties'
    return None

def test_dual_gauss_1d_nounc():
    # Generate data from two Gaussians, recover mean and variance
    ndata= 3001
    amp_true= 0.3
    assign= numpy.random.binomial(1,1.-amp_true,ndata)
    ydata= numpy.zeros((ndata,1))
    ydata[assign==0,0]= numpy.random.normal(size=numpy.sum(assign==0))-2.
    ydata[assign==1,0]= numpy.random.normal(size=numpy.sum(assign==1))*2.+1.
    ycovar= numpy.zeros_like(ydata)
    # initialize fit
    K= 2
    initamp= numpy.ones(K)/float(K)
    initmean= numpy.zeros((K,1))
    initcovar= numpy.zeros((K,1,1))
    for kk in range(K):
        initmean[kk]= numpy.mean(ydata)+numpy.random.normal()*numpy.std(ydata)
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 6./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover correct amp for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover correct mean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover correct amp for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover correct mean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties'
    return None

def test_dual_gauss_1d_constunc():
    # Generate data from two Gaussians, recover mean and variance
    ndata= 3001
    amp_true= 0.3
    assign= numpy.random.binomial(1,1.-amp_true,ndata)
    ydata= numpy.zeros((ndata,1))
    ydata[assign==0,0]= numpy.random.normal(size=numpy.sum(assign==0))-2.
    ydata[assign==1,0]= numpy.random.normal(size=numpy.sum(assign==1))*2.+1.
    ycovar= numpy.ones_like(ydata)*0.25
    ydata+= numpy.atleast_2d(numpy.random.normal(size=ndata)).T\
        *numpy.sqrt(ycovar)
    # initialize fit
    K= 2
    initamp= numpy.ones(K)/float(K)
    initmean= numpy.zeros((K,1))
    initcovar= numpy.zeros((K,1,1))
    numpy.random.uniform() # for the current rng, initmean would be very similar without this...
    for kk in range(K):
        initmean[kk]= numpy.mean(ydata)+numpy.random.normal()*numpy.std(ydata)
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 10./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover correct amp for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover correct mean for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/ constant uncertainties'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover correct amp for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover correct mean for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/ constant uncertainties'
    return None

def test_dual_gauss_1d_varunc():
    # Generate data from two Gaussians, recover mean and variance
    ndata= 3001
    amp_true= 0.3
    assign= numpy.random.binomial(1,1.-amp_true,ndata)
    ydata= numpy.zeros((ndata,1))
    ydata[assign==0,0]= numpy.random.normal(size=numpy.sum(assign==0))-2.
    ydata[assign==1,0]= numpy.random.normal(size=numpy.sum(assign==1))*2.+1.
    ycovar= numpy.ones_like(ydata)*\
        numpy.atleast_2d(numpy.random.uniform(size=ndata)).T
    ydata+= numpy.atleast_2d(numpy.random.normal(size=ndata)).T\
        *numpy.sqrt(ycovar)
    # initialize fit
    K= 2
    initamp= numpy.ones(K)/float(K)
    initmean= numpy.zeros((K,1))
    initcovar= numpy.zeros((K,1,1))
    numpy.random.uniform() # for the current rng, initmean would be very similar without this...
    for kk in range(K):
        initmean[kk]= numpy.mean(ydata)+numpy.random.normal()*numpy.std(ydata)
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 10./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover correct amp for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover correct mean for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/ constant uncertainties'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover correct amp for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover correct mean for dual Gaussian w/ constant uncertainties'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/ constant uncertainties'
    return None

