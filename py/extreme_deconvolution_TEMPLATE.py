import os, os.path
import ctypes
import ctypes.util
import numpy as nu
from numpy.ctypeslib import ndpointer
#Find and load the library
_lib = None
_libname = ctypes.util.find_library('libextremedeconvolution.so')
if _libname:
    _lib = ctypes.CDLL(_libname)
if _lib is None: #Hack
    p = os.path.join(TEMPLATE_LIBRARY_PATH,'libextremedeconvolution.so')
    if os.path.exists(p):
        _lib = ctypes.CDLL(p)
if _lib is None:
        raise IOError('libextremedeconvolution.so library not found')

def _fix2chararray(fix,ngauss):
    """Internal function to process the fix* inputs"""
    if fix == None:
        fix= [chr(False) for kk in range(ngauss)]
    else: #fix is set
        try:
            if len(fix) == 1 and ngauss != 1:
                fixamp= [chr(fixamp[0]) for kk in range(ngauss)]
            else: #We assume all values are set
                fix= [chr(fix[kk]) for kk in range(ngauss)]
        except TypeError: #fixamp == Bool
            fix= [chr(fix) for kk in range(ngauss)]
    return fix

def extreme_deconvolution(ydata,ycovar,
                          xamp,xmean,xcovar,
                          projection=None,
                          fixamp=None,fixmean=None,fixcovar=None,
                          tol=1.e-6,maxiter=long(1e9),w=0.,logfile=None,
                          splitnmerge=0,maxsnm=False,likeonly=False):
    """
    NAME:
       extreme_deconvolution
    PURPOSE:
       run the underlying C-extreme-deconvolution code
    INPUT:
       ydata - [ndata,dy] numpy array of observed quantities
       ycovar - [ndata,dy,dy] numpy array of observational error covariances
       xamp - [ngauss] numpy array of initial amplitudes
       xmean - [ngauss,dx] numpy array of initial means
       xcovar - [ngauss,dx,dx] numpy array of initial covariances
    OPTIONAL INPUTS:
       projection - [ndata,dy,dx] numpy array of projection matrices
       fixamp - (default=None) None, True/False, or list of bools
       fixmean - (default=None) None, True/False, or list of bools
       fixcovar - (default=None) None, True/False, or list of bools
       tol - (double, default=1.e-6) tolerance for convergence
       maxiter - (long, default= 10**9) maximum number of iterations to perform
       w - (double, default=0.) covariance regularization parameter
            (of the conjugate prior)
       logfile - basename for several logfiles (_c.log has output from
                 the c-routine; _loglike.log has the log likelihood path of
                 all the accepted routes, i.e. only parts which increase
                 the likelihood are included, during splitnmerge)
       splitnmerge - (int, default=0) depth to go down the splitnmerge path
       maxsnm - (Bool, default=False) use the maximum number of split 'n'
                 merge steps, K*(K-1)*(K-2)/2
       likeonly - (Bool, default=False) only compute the total log
                   likelihood of the data
    OUTPUT:
       avgloglikedata after convergence,
       +updated xamp, xmean, xcovar
    HISTORY:
       2010-02-10 - Written - Bovy (NYU)
    DOCTEST:
    >>> import numpy as nu
    >>> nu.random.seed(seed=-1)
    >>> ndata=200
    >>> ydata=nu.reshape(nu.random.normal(1,1,ndata),(ndata,1))
    >>> ycovar= nu.reshape(nu.ones(ndata)*.01,(ndata,1,1))
    >>> projection= nu.zeros((ndata,1,2))
    >>> xamp= nu.ones(2)/2.
    >>> xmean= nu.random.random((2,2))
    >>> xcovar= nu.random.random((2,2,2))/10.
    >>> for ii in range(ndata):
    ...     if (ii % 2) == 0:
    ...             projection[ii,:,:]= nu.array([1.,0.])
    ...     else:
    ...             projection[ii,:,:]= nu.array([0.,1.])
    ...
    >>> print xmean
    [[ 0.6613586   0.75460518]
     [ 0.59746367  0.17107033]]
    >>> extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,projection=projection)
    -1.3594338204401231
    >>> print xmean
    [[ 0.0975456   0.03640963]
     [ 1.41598304  1.15101783]]
     """
    ndata= ydata.shape[0]
    dataDim= ydata.shape[1]
    ngauss= len(xamp)
    gaussDim= xmean.shape[1]

    fixamp= _fix2chararray(fixamp,ngauss)
    fixmean= _fix2chararray(fixmean,ngauss)
    fixcovar= _fix2chararray(fixcovar,ngauss)

    avgloglikedata= ctypes.pointer(ctypes.c_double(0.))

    if logfile == None:
        clog= ''
        clog2= ''
        n_clog= 0
        n_clog2= 0
    else:
        clog= logfile + '_c.log'
        n_clog= len(clog)
        clog2= logfile + '_loglike.log'
        n_clog2= len(clog2)

    if maxsnm:
        splitnmerge = long(ngauss*(ngauss-1)*(ngauss-2)/2)

    if projection == None:
        noprojection= True
        projection= nu.zeros(1)
    else:
        noprojection= False
        
    ndarrayFlags= ('C_CONTIGUOUS','WRITEABLE')
    exdeconvFunc= _lib.proj_gauss_mixtures_IDL
    exdeconvFunc.argtypes = [ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ctypes.c_int,
                             ctypes.c_int,
                             ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ctypes.c_int,
                             ctypes.c_int,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.POINTER(ctypes.c_double),
                             ctypes.c_double,
                             ctypes.c_int,
                             ctypes.c_char,
                             ctypes.c_double,
                             ctypes.c_char_p,
                             ctypes.c_int,
                             ctypes.c_int,
                             ctypes.c_char_p,
                             ctypes.c_int,
                             ctypes.c_char]
                                             
                                             
    exdeconvFunc(ydata,
                 ycovar,
                 projection,
                 ctypes.c_int(ndata),
                 ctypes.c_int(dataDim),
                 xamp,
                 xmean,
                 xcovar,
                 ctypes.c_int(gaussDim),
                 ctypes.c_int(ngauss),
                 ctypes.byref(ctypes.c_char(fixamp[0])),
                 ctypes.byref(ctypes.c_char(fixmean[0])),
                 ctypes.byref(ctypes.c_char(fixcovar[0])),
                 avgloglikedata,
                 ctypes.c_double(tol),
                 ctypes.c_int(maxiter),
                 ctypes.c_char(chr(likeonly)),
                 ctypes.c_double(w),
                 ctypes.create_string_buffer(clog),
                 ctypes.c_int(n_clog),
                 ctypes.c_int(splitnmerge),
                 ctypes.create_string_buffer(clog2),
                 ctypes.c_int(n_clog2),
                 ctypes.c_char(chr(noprojection)))
    return avgloglikedata.contents.value

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
