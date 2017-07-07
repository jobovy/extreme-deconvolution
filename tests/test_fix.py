# test_fix.py: test that fixing parameters works as expected
import numpy
numpy.random.seed(2)
from extreme_deconvolution import extreme_deconvolution

# Very basic tests in 1D with K=1, cannot do amp
def test_fixmean_single_gauss_1d_nounc():
    # Generate data from a single Gaussian, recover variance, fixing mean
    ndata= 3001
    ydata= numpy.atleast_2d(numpy.random.normal(size=ndata)).T+1.
    ycovar= numpy.zeros_like(ydata)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.array([[1.]])
    initcovar= numpy.atleast_3d(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixmean=True)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean-1.) < 10.**-10., 'XD did not fixmean for single Gaussian'
    assert numpy.fabs(initcovar-1.) < tol, 'XD does not recover correct variance for single Gaussian w/o uncertainties, fixing mean'
    return None

def test_fixcovar_single_gauss_1d_nounc():
    # Generate data from a single Gaussian, recover variance, fixing mean
    ndata= 3001
    ydata= numpy.atleast_2d(numpy.random.normal(size=ndata)).T+1.
    ycovar= numpy.zeros_like(ydata)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata)+numpy.std(ydata))
    initcovar= numpy.array([[[1.]]])
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixcovar=True)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean-1.) < tol, 'XD does not recover mean for single Gaussian, fixing covar'
    assert numpy.fabs(initcovar-1.) < 10.**-10., 'XD did not fixcovar for single Gaussian w/o uncertainties'
    return None

# More complicated cases for K=2
def test_fixamp_dual_gauss_1d_nounc():
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
    initamp= numpy.array([amp_true,1.-amp_true])
    initmean= numpy.array([[-1.],[2.]])
    initcovar= numpy.zeros((K,1,1))
    for kk in range(K):
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixamp=True)
    # Test
    tol= 12./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < 10.**-10., 'XD did not fixamp for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover correct mean for dual Gaussian w/o uncertainties, fixing amp'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing amp'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < 10.**-10., 'XD did not fixamp for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover correct mean for dual Gaussian w/o uncertainties, fixing amp'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing amp'
    return None

def test_fixamp_alt_dual_gauss_1d_nounc():
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
    initamp= numpy.array([amp_true,1.-amp_true])
    initmean= numpy.array([[-1.],[2.]])
    initcovar= numpy.zeros((K,1,1))
    numpy.random.uniform() # hack to get diff init
    for kk in range(K):
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixamp=[True,False]) # should be same as =True
    # Test
    tol= 12./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < 10.**-10., 'XD did not fixamp for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover correct mean for dual Gaussian w/o uncertainties, fixing amp'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing amp'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < 10.**-10., 'XD did not fixamp for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover correct mean for dual Gaussian w/o uncertainties, fixing amp'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing amp'
    return None

def test_fixmean_dual_gauss_1d_nounc():
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
        if kk == 0:
            initmean[kk]= -2.
        elif kk == 1:
            initmean[kk]= 1.
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixmean=True)
    # Test
    tol= 12./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing mean'
    assert numpy.fabs(initmean[first]--2.) < 10.**-10., 'XD did not fixmean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing mean'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing mean'
    assert numpy.fabs(initmean[second]-1.) < 10.**-10., 'XD did not fixmean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing covar'
    return None

def test_fixmean_fixone_dual_gauss_1d_nounc():
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
        if kk == 0:
            initmean[kk]= -2.
        elif kk == 1:
            initmean[kk]= 1.5
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixmean=[True,False])
    # Test
    tol= 12./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing one mean'
    assert numpy.fabs(initmean[first]--2.) < 10.**-10., 'XD did not fixmean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing one mean'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing one mean'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover mean  for dual Gaussian w/o uncertainties, fixing one mean'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing one mean'
    return None

def test_fixcovar_dual_gauss_1d_nounc():
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
    initmean= numpy.array([[-1.],[2.]])
    initcovar= numpy.zeros((K,1,1))
    for kk in range(K):
        if kk == 0:
            initcovar[kk]= 1.
        elif kk == 1:
            initcovar[kk]= 4.
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixcovar=True)
    # Test
    tol= 12./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing covar'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover mean for dual Gaussian w/o uncertainties, fixing covar'
    assert numpy.fabs(initcovar[first]-1.) < tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing mean'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing covar'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover mean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[second]-4.) < 10.**-10., 'XD did not fixcovar for dual Gaussian w/o uncertainties'
    return None

def test_fixcovar_fixone_dual_gauss_1d_nounc():
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
    initmean= numpy.array([[-1.],[2.]])
    initcovar= numpy.zeros((K,1,1))
    for kk in range(K):
        if kk == 0:
            initcovar[kk]= 1.
        elif kk == 1:
            initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          fixcovar=[True,False])
    # Test
    tol= 12./numpy.sqrt(ndata)
    first= initamp < 0.5
    assert numpy.fabs(initamp[first]-amp_true) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing one covar'
    assert numpy.fabs(initmean[first]--2.) < tol, 'XD does not recover mean for dual Gaussian w/o uncertainties'
    assert numpy.fabs(initcovar[first]-1.) < 10.**-10., 'XD did not fixcovar for dual Gaussian w/o uncertainties'
    second= initamp >= 0.5
    assert numpy.fabs(initamp[second]-(1.-amp_true)) < tol, 'XD does not recover amp for dual Gaussian w/o uncertainties, fixing one covar'
    assert numpy.fabs(initmean[second]-1.) < 2.*tol, 'XD does not recover mean  for dual Gaussian w/o uncertainties, fixing one covar'
    assert numpy.fabs(initcovar[second]-4.) < 2.*tol, 'XD does not recover correct variance for dual Gaussian w/o uncertainties, fixing one covar'
    return None

