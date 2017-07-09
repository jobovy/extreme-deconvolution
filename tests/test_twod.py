# test_twod.py: tests of XD in 2D
import numpy
numpy.random.seed(2)
from extreme_deconvolution import extreme_deconvolution

def test_single_gauss_2d_nounc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    ydata= numpy.random.normal(size=(ndata,2))+numpy.array([[1.,2.]])
    ycovar= numpy.zeros_like(ydata)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata,axis=0)\
                                   +numpy.std(ydata,axis=0))
    initcovar= numpy.atleast_3d(numpy.cov(ydata,rowvar=False)).T
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean[0,0]-1.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initmean[0,1]-2.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,0]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,1,1]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,1]-0.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    return None

def test_single_gauss_2d_diagunc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    ydata= numpy.random.normal(size=(ndata,2))+numpy.array([[1.,2.]])
    ycovar= numpy.ones_like(ydata)\
        *numpy.random.uniform(size=(ndata,2))/2.
    ydata+= numpy.random.normal(size=(ndata,2))*numpy.sqrt(ycovar)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata,axis=0)\
                                   +numpy.std(ydata,axis=0))
    initcovar= numpy.atleast_3d(numpy.cov(ydata,rowvar=False)).T
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean[0,0]-1.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initmean[0,1]-2.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,0]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,1,1]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,1]-0.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    return None

def test_single_gauss_2d_offdiagunc():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    ydata= numpy.random.normal(size=(ndata,2))+numpy.array([[1.,2.]])
    tycovar= numpy.ones_like(ydata)\
        *numpy.random.uniform(size=(ndata,2))/2.
    ydata+= numpy.random.normal(size=(ndata,2))*numpy.sqrt(tycovar)
    ycovar= numpy.empty((ndata,2,2))
    for ii in range(ndata):
        ycovar[ii]= numpy.diag(tycovar[ii])
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata,axis=0)\
                                   +numpy.std(ydata,axis=0))
    initcovar= numpy.atleast_3d(numpy.cov(ydata,rowvar=False)).T
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean[0,0]-1.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initmean[0,1]-2.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,0]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,1,1]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,1]-0.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    return None

def test_single_gauss_2d_diagunc_proj():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    tydata= numpy.random.normal(size=(ndata,2))+numpy.array([[1.,2.]])
    # Randomly project
    ydata= numpy.zeros((ndata,1))
    proj_x= numpy.random.binomial(2,0.5,ndata)
    ydata[proj_x==0,0]= tydata[proj_x==0,0]
    ydata[proj_x==1,0]= tydata[proj_x==1,1]
    ydata[proj_x==2,0]= tydata[proj_x==2,0]+tydata[proj_x==2,1]
    projection= numpy.empty((ndata,1,2))
    projection[proj_x==0]= numpy.array([[[1.,0.]]])
    projection[proj_x==1]= numpy.array([[[0.,1.]]])
    projection[proj_x==2]= numpy.array([[[1.,1.]]])
    ycovar= numpy.ones_like(ydata)*\
        numpy.atleast_2d(numpy.random.uniform(size=ndata)).T
    ydata+= numpy.atleast_2d(numpy.random.normal(size=ndata)).T\
        *numpy.sqrt(ycovar)
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.array([[0.,1.]])
    initcovar= numpy.array([[[2.,-1.],[-1.,3.]]])
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          projection=projection)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean[0,0]-1.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initmean[0,1]-2.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,0]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,1,1]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,1]-0.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    return None

def test_single_gauss_2d_offdiagunc_proj():
    # Generate data from a single Gaussian, recover mean and variance
    ndata= 3001
    ydata= numpy.random.normal(size=(ndata,2))+numpy.array([[1.,2.]])
    # For half of the points, x -> x+y
    proj= numpy.random.uniform(size=ndata) > 0.5
    ydata[proj,0]= ydata[proj,0]+ydata[proj,1]
    projection= numpy.empty((ndata,2,2))
    projection[proj]= numpy.array([[1.,1.],[0.,1.]])
    projection[True^proj]= numpy.array([[1.,0.],[0.,1.]])
    tycovar= numpy.ones_like(ydata)\
        *numpy.random.uniform(size=(ndata,2))/2.
    ydata+= numpy.random.normal(size=(ndata,2))*numpy.sqrt(tycovar)
    ycovar= numpy.empty((ndata,2,2))
    for ii in range(ndata):
        ycovar[ii]= numpy.diag(tycovar[ii])
    # initialize fit
    K= 1
    initamp= numpy.ones(K)
    initmean= numpy.atleast_2d(numpy.mean(ydata,axis=0)\
                                   +numpy.std(ydata,axis=0))
    initcovar= numpy.atleast_3d(numpy.cov(ydata,rowvar=False)).T
    # Run XD
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          projection=projection)
    # Test
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean[0,0]-1.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initmean[0,1]-2.) < tol, 'XD does not recover correct mean for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,0]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,1,1]-1.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    assert numpy.fabs(initcovar[0,0,1]-0.) < tol, 'XD does not recover correct variance for single Gaussian in 2D w/o uncertainties'
    return None

