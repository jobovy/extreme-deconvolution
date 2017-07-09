# test_log.py: Test that things are getting logged
import os, os.path
import numpy
numpy.random.seed(2)
from extreme_deconvolution import extreme_deconvolution

def test_single_gauss_1d_varunc_log():
    # Same as in test_oned, but now also log
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
    logfile= 'test_log'
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          logfile=logfile)
    # First test that fit worked
    tol= 10./numpy.sqrt(ndata)
    assert numpy.fabs(initmean-0.) < tol, 'XD does not recover correct mean for single Gaussian w/ uncertainties'
    assert numpy.fabs(initcovar-1.) < tol, 'XD does not recover correct variance for single Gaussian w/ uncertainties'
    # Now test that the logfiles exist
    assert os.path.exists(logfile+'_c.log'), 'XD did not produce _c.log logfile when asked'
    num_lines= sum(1 for line in open(logfile+'_c.log'))
    assert num_lines > 0, "XD logfile _c.log appears to be empty, but shouldn't be"
    assert os.path.exists(logfile+'_loglike.log'), 'XD did not produce _loglike.log logfile when asked'
    num_lines= sum(1 for line in open(logfile+'_loglike.log'))
    assert num_lines > 0, "XD logfile _loglike.log appears to be empty, but shouldn't be"
    os.remove(logfile+'_c.log')
    os.remove(logfile+'_loglike.log')
    return None

def test_triple_gauss_1d_varunc_snm_log():
    # Like in oned, but also log
    ndata= 3001
    amp_true= [0.3,0.1,0.6]
    assign= numpy.random.choice(numpy.arange(3),p=amp_true,size=ndata)
    ydata= numpy.zeros((ndata,1))
    ydata[assign==0,0]= numpy.random.normal(size=numpy.sum(assign==0))-4.
    ydata[assign==1,0]= numpy.random.normal(size=numpy.sum(assign==1))*2.+1.
    ydata[assign==2,0]= numpy.random.normal(size=numpy.sum(assign==2))*1.5+8.
    ycovar= numpy.ones_like(ydata)*\
        numpy.atleast_2d(numpy.random.uniform(size=ndata)).T
    ydata+= numpy.atleast_2d(numpy.random.normal(size=ndata)).T\
        *numpy.sqrt(ycovar)
    # initialize fit
    K= 3
    initamp= numpy.ones(K)/float(K)
    initmean= numpy.array([[-1.],[0.],[1.]])
    initcovar= numpy.zeros((K,1,1))
    for kk in range(K):
        initcovar[kk]= numpy.mean(3.*numpy.var(ydata))
    # Run XD
    logfile= 'test_log'
    extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar,
                          maxsnm=True,logfile=logfile)
    # Test
    tol= 25./numpy.sqrt(ndata)
    first= initamp > 0.5
    assert numpy.fabs(initamp[first]-amp_true[2]) < tol, 'XD does not recover correct amp for triple Gaussian w/ uncertainties'
    assert numpy.fabs(initmean[first]-8.) < tol, 'XD does not recover correct mean for triple Gaussian w/ uncertainties'
    assert numpy.fabs(initcovar[first]-1.5**2.) < tol, 'XD does not recover correct variance for triple Gaussian w/ uncertainties'
    second= (initamp <= 0.5)*(initamp > 0.2)
    assert numpy.fabs(initamp[second]-amp_true[0]) < tol, 'XD does not recover correct amp for triple Gaussian w/  uncertainties'
    assert numpy.fabs(initmean[second]--4.) < 2.*tol, 'XD does not recover correct mean for triple Gaussian w/  uncertainties'
    assert numpy.fabs(initcovar[second]-1.) < 2.*tol, 'XD does not recover correct variance for triple Gaussian w/  uncertainties'
    third= (initamp <= 0.2)
    assert numpy.fabs(initamp[third]-amp_true[1]) < tol, 'XD does not recover correct amp for triple Gaussian w/  uncertainties'
    assert numpy.fabs(initmean[third]-1.) < 4.*tol, 'XD does not recover correct mean for triple Gaussian w/  uncertainties'
    assert numpy.fabs(initcovar[third]-4.) < 6.*tol, 'XD does not recover correct variance for triple Gaussian w/  uncertainties'
    # Now test that the logfiles exist
    assert os.path.exists(logfile+'_c.log'), 'XD did not produce _c.log logfile when asked'
    num_lines= sum(1 for line in open(logfile+'_c.log'))
    assert num_lines > 0, "XD logfile _c.log appears to be empty, but shouldn't be"
    assert os.path.exists(logfile+'_loglike.log'), 'XD did not produce _loglike.log logfile when asked'
    num_lines= sum(1 for line in open(logfile+'_loglike.log'))
    assert num_lines > 0, "XD logfile _loglike.log appears to be empty, but shouldn't be"
    os.remove(logfile+'_c.log')
    os.remove(logfile+'_loglike.log')
    return None

