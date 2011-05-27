import os, os.path, platform
import ctypes
import ctypes.util
import numpy as nu
from numpy.ctypeslib import ndpointer
#Find and load the library
_lib = None
if platform.system()=='Darwin':
    _libraryname= 'libextremedeconvolution.dylib'
else:
    _libraryname= 'libextremedeconvolution.so'
_libname = ctypes.util.find_library(_libraryname)
if _libname:
    _lib = ctypes.CDLL(_libname)
if _lib is None: #Hack
    p = os.path.join('/usr/local/lib/',_libraryname)
    if os.path.exists(p):
        _lib = ctypes.CDLL(p)
if _lib is None:
        raise IOError(_libraryname+' library not found')

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
                          weight=None,
                          fixamp=None,fixmean=None,fixcovar=None,
                          tol=1.e-6,maxiter=long(1e9),w=0.,logfile=None,
                          splitnmerge=0,maxsnm=False,likeonly=False,
                          logweight=False):
    """
    NAME:
       extreme_deconvolution
    PURPOSE:
       run the underlying C-extreme-deconvolution code
    INPUT:
       ydata - [ndata,dy] numpy array of observed quantities
       ycovar - [ndata,dy(,dy)] numpy array of observational error covariances
                (if [ndata,dy] then the error correlations are assumed to vanish)
       xamp - [ngauss] numpy array of initial amplitudes
       xmean - [ngauss,dx] numpy array of initial means
       xcovar - [ngauss,dx,dx] numpy array of initial covariances
    OPTIONAL INPUTS:
       projection - [ndata,dy,dx] numpy array of projection matrices
       weight - [ndata] numpy array of weights to be applied to the data points
       logweight - (bool, default=False) if True, weight is actually
                   log(weight)
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
    >>> nu.random.seed(seed=1)
    >>> ydata= nu.array([[  2.62434536e+00],
    ...                  [  3.88243586e-01],
    ...                  [  4.71828248e-01],
    ...                  [ -7.29686222e-02],
    ...                  [  1.86540763e+00],
    ...                  [ -1.30153870e+00],
    ...                  [  2.74481176e+00],
    ...                  [  2.38793099e-01],
    ...                  [  1.31903910e+00],
    ...                  [  7.50629625e-01],
    ...                  [  2.46210794e+00],
    ...                  [ -1.06014071e+00],
    ...                  [  6.77582796e-01],
    ...                  [  6.15945645e-01],
    ...                  [  2.13376944e+00],
    ...                  [ -9.98912673e-02],
    ...                  [  8.27571792e-01],
    ...                  [  1.22141582e-01],
    ...                  [  1.04221375e+00],
    ...                  [  1.58281521e+00],
    ...                  [ -1.00619177e-01],
    ...                  [  2.14472371e+00],
    ...                  [  1.90159072e+00],
    ...                  [  1.50249434e+00],
    ...                  [  1.90085595e+00],
    ...                  [  3.16272141e-01],
    ...                  [  8.77109774e-01],
    ...                  [  6.42305657e-02],
    ...                  [  7.32111920e-01],
    ...                  [  1.53035547e+00],
    ...                  [  3.08339248e-01],
    ...                  [  6.03246473e-01],
    ...                  [  3.12827300e-01],
    ...                  [  1.54794359e-01],
    ...                  [  3.28753869e-01],
    ...                  [  9.87335401e-01],
    ...                  [ -1.17310349e-01],
    ...                  [  1.23441570e+00],
    ...                  [  2.65980218e+00],
    ...                  [  1.74204416e+00],
    ...                  [  8.08164448e-01],
    ...                  [  1.12371036e-01],
    ...                  [  2.52841706e-01],
    ...                  [  2.69245460e+00],
    ...                  [  1.05080775e+00],
    ...                  [  3.63004353e-01],
    ...                  [  1.19091548e+00],
    ...                  [  3.10025514e+00],
    ...                  [  1.12015895e+00],
    ...                  [  1.61720311e+00],
    ...                  [  1.30017032e+00],
    ...                  [  6.47750154e-01],
    ...                  [ -1.42518198e-01],
    ...                  [  6.50657278e-01],
    ...                  [  7.91105767e-01],
    ...                  [  1.58662319e+00],
    ...                  [  1.83898341e+00],
    ...                  [  1.93110208e+00],
    ...                  [  1.28558733e+00],
    ...                  [  1.88514116e+00],
    ...                  [  2.45602059e-01],
    ...                  [  2.25286816e+00],
    ...                  [  1.51292982e+00],
    ...                  [  7.01907165e-01],
    ...                  [  1.48851815e+00],
    ...                  [  9.24428287e-01],
    ...                  [  2.13162939e+00],
    ...                  [  2.51981682e+00],
    ...                  [  3.18557541e+00],
    ...                  [ -3.96496335e-01],
    ...                  [ -4.44113805e-01],
    ...                  [  4.95534137e-01],
    ...                  [  1.16003707e+00],
    ...                  [  1.87616892e+00],
    ...                  [  1.31563495e+00],
    ...                  [ -1.02220122e+00],
    ...                  [  6.93795987e-01],
    ...                  [  1.82797464e+00],
    ...                  [  1.23009474e+00],
    ...                  [  1.76201118e+00],
    ...                  [  7.77671857e-01],
    ...                  [  7.99241931e-01],
    ...                  [  1.18656139e+00],
    ...                  [  1.41005165e+00],
    ...                  [  1.19829972e+00],
    ...                  [  1.11900865e+00],
    ...                  [  3.29337714e-01],
    ...                  [  1.37756379e+00],
    ...                  [  1.12182127e+00],
    ...                  [  2.12948391e+00],
    ...                  [  2.19891788e+00],
    ...                  [  1.18515642e+00],
    ...                  [  6.24715050e-01],
    ...                  [  3.61269593e-01],
    ...                  [  1.42349435e+00],
    ...                  [  1.07734007e+00],
    ...                  [  6.56146324e-01],
    ...                  [  1.04359686e+00],
    ...                  [  3.79999156e-01],
    ...                  [  1.69803203e+00],
    ...                  [  5.52871435e-01],
    ...                  [  2.22450770e+00],
    ...                  [  1.40349164e+00],
    ...                  [  1.59357852e+00],
    ...                  [ -9.49118457e-02],
    ...                  [  1.16938243e+00],
    ...                  [  1.74055645e+00],
    ...                  [  4.62993982e-02],
    ...                  [  7.33781494e-01],
    ...                  [  1.03261455e+00],
    ...                  [ -3.73117320e-01],
    ...                  [  1.31515939e+00],
    ...                  [  1.84616065e+00],
    ...                  [  1.40484059e-01],
    ...                  [  1.35054598e+00],
    ...                  [ -3.12283411e-01],
    ...                  [  9.61304491e-01],
    ...                  [ -6.15772355e-01],
    ...                  [  2.12141771e+00],
    ...                  [  1.40890054e+00],
    ...                  [  9.75383044e-01],
    ...                  [  2.24838381e-01],
    ...                  [  2.27375593e+00],
    ...                  [  2.96710175e+00],
    ...                  [ -8.57981864e-01],
    ...                  [  2.23616403e+00],
    ...                  [  2.62765075e+00],
    ...                  [  1.33801170e+00],
    ...                  [ -1.99268032e-01],
    ...                  [  1.86334532e+00],
    ...                  [  8.19079698e-01],
    ...                  [  3.96079372e-01],
    ...                  [ -2.30058136e-01],
    ...                  [  1.55053750e+00],
    ...                  [  1.79280687e+00],
    ...                  [  3.76469270e-01],
    ...                  [  1.52057634e+00],
    ...                  [ -1.44341390e-01],
    ...                  [  1.80186103e+00],
    ...                  [  1.04656730e+00],
    ...                  [  8.13430228e-01],
    ...                  [  8.98254127e-01],
    ...                  [  1.86888616e+00],
    ...                  [  1.75041164e+00],
    ...                  [  1.52946532e+00],
    ...                  [  1.13770121e+00],
    ...                  [  1.07782113e+00],
    ...                  [  1.61838026e+00],
    ...                  [  1.23249456e+00],
    ...                  [  1.68255141e+00],
    ...                  [  6.89883226e-01],
    ...                  [ -1.43483776e+00],
    ...                  [  2.03882460e+00],
    ...                  [  3.18697965e+00],
    ...                  [  1.44136444e+00],
    ...                  [  8.99844767e-01],
    ...                  [  8.63555256e-01],
    ...                  [  8.80945812e-01],
    ...                  [  1.01740941e+00],
    ...                  [ -1.22018729e-01],
    ...                  [  4.82905542e-01],
    ...                  [  2.97317235e-03],
    ...                  [  1.24879916e+00],
    ...                  [  7.03358848e-01],
    ...                  [  1.49521132e+00],
    ...                  [  8.25296840e-01],
    ...                  [  1.98633519e+00],
    ...                  [  1.21353390e+00],
    ...                  [  3.19069973e+00],
    ...                  [ -8.96360923e-01],
    ...                  [  3.53083312e-01],
    ...                  [  1.90148689e+00],
    ...                  [  3.52832571e+00],
    ...                  [  7.51365222e-01],
    ...                  [  1.04366899e+00],
    ...                  [  7.73685757e-01],
    ...                  [  2.33145711e+00],
    ...                  [  7.12692137e-01],
    ...                  [  1.68006984e+00],
    ...                  [  6.80198401e-01],
    ...                  [ -2.72558755e-01],
    ...                  [  1.31354772e+00],
    ...                  [  1.50318481e+00],
    ...                  [  2.29322588e+00],
    ...                  [  8.89552974e-01],
    ...                  [  3.82637936e-01],
    ...                  [  1.56276110e+00],
    ...                  [  1.24073709e+00],
    ...                  [  1.28066508e+00],
    ...                  [  9.26887296e-01],
    ...                  [  2.16033857e+00],
    ...                  [  1.36949272e+00],
    ...                  [  2.90465871e+00],
    ...                  [  2.11105670e+00],
    ...                  [  1.65904980e+00],
    ...                  [ -6.27438341e-01],
    ...                  [  1.60231928e+00],
    ...                  [  1.42028220e+00],
    ...                  [  1.81095167e+00],
    ...                  [  2.04444209e+00]])
    >>> ndata= len(ydata)
    >>> ycovar= nu.reshape(nu.ones(ndata)*.01,(ndata,1))
    >>> projection= nu.zeros((ndata,1,2))
    >>> xamp= nu.ones(2)/2.
    >>> xmean= nu.array([[ 0.86447943,0.322681  ],
    ...                  [ 0.67078879,0.45087394]])
    >>> xcovar= nu.array([[[ 0.03821028, 0.04108113],
    ...                   [ 0.04014796,  0.03173839]],
    ...                   [[ 0.06219194,  0.04302473],
    ...                   [ 0.09738021,  0.06778009]]])
    >>> for ii in range(ndata):
    ...     if (ii % 2) == 0:
    ...             projection[ii,:,:]= nu.array([1.,0.])
    ...     else:
    ...             projection[ii,:,:]= nu.array([0.,1.])
    ...
    >>> l= extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,projection=projection)
    >>> assert (l--1.3114744655258121)**2. < 10.**-8
    >>> assert (xmean[0,0]-2.30368235)**2. < 10.**-5
    >>> assert (xmean[0,1]-1.70701517)**2. < 10.**-5
    >>> assert (xmean[1,0]-1.08009397)**2. < 10.**-5
    >>> assert (xmean[1,1]-0.8888667)**2. < 10.**-5
    >>> ydata= nu.asfortranarray(ydata)
    >>> ydata.flags['C_CONTIGUOUS']
    False
    >>> ydata.flags['F_CONTIGUOUS']
    True
    >>> l= extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,projection=projection)
    >>> assert (l--1.3114744655258121)**2. < 10.**-8
    >>> assert (xmean[0,0]-2.30368235)**2. < 10.**-5
    >>> assert (xmean[0,1]-1.70701517)**2. < 10.**-5
    >>> assert (xmean[1,0]-1.08009397)**2. < 10.**-5
    >>> assert (xmean[1,1]-0.8888667)**2. < 10.**-5
    >>> ydata.flags['C_CONTIGUOUS']
    False
    >>> ydata.flags['F_CONTIGUOUS']
    True
    """
    ndata= ydata.shape[0]
    dataDim= ydata.shape[1]
    ngauss= len(xamp)
    gaussDim= xmean.shape[1]

    if len(ycovar.shape) == 2:
        diagerrors= True
    else:
        diagerrors= False

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
        
    if weight == None:
        noweight= True
        logweights= nu.zeros(1)
    elif not logweight:
        noweight= False
        logweights= nu.log(weight)
    else:
        noweight= False
        logweights= weight
        
    ndarrayFlags= ('C_CONTIGUOUS','WRITEABLE')
    exdeconvFunc= _lib.proj_gauss_mixtures_IDL
    exdeconvFunc.argtypes = [ndpointer(dtype=nu.float64,flags=ndarrayFlags),
                             ndpointer(dtype=nu.float64,flags=ndarrayFlags),
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
                             ctypes.c_char,
                             ctypes.c_char,
                             ctypes.c_char]
                                             
    #Requirements, first store old order
    f_cont= [ydata.flags['F_CONTIGUOUS'],
             ycovar.flags['F_CONTIGUOUS'],
             projection.flags['F_CONTIGUOUS'],
             logweights.flags['F_CONTIGUOUS'],
             xamp.flags['F_CONTIGUOUS'],
             xmean.flags['F_CONTIGUOUS'],
             xcovar.flags['F_CONTIGUOUS']]
    ydata= nu.require(ydata,dtype=nu.float64,requirements=['C','W'])
    ycovar= nu.require(ycovar,dtype=nu.float64,requirements=['C','W'])
    projection= nu.require(projection,dtype=nu.float64,requirements=['C','W'])
    logweights= nu.require(logweights,dtype=nu.float64,requirements=['C','W'])
    xamp= nu.require(xamp,dtype=nu.float64,requirements=['C','W'])
    xmean= nu.require(xmean,dtype=nu.float64,requirements=['C','W'])
    xcovar= nu.require(xcovar,dtype=nu.float64,requirements=['C','W'])

    exdeconvFunc(ydata,
                 ycovar,
                 projection,
                 logweights,
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
                 ctypes.c_char(chr(noprojection)),
                 ctypes.c_char(chr(diagerrors)),
                 ctypes.c_char(chr(noweight)))
    #Reset input arrays
    if f_cont[0]: ydata= nu.asfortranarray(ydata)
    if f_cont[1]: ycovar= nu.asfortranarray(ycovar)
    if f_cont[2]: projection= nu.asfortranarray(projection)
    if f_cont[3]: logweights= nu.asfortranarray(logweights)
    if f_cont[4]: xamp=  nu.asfortranarray(xamp)
    if f_cont[5]: xmean=  nu.asfortranarray(xmean)
    if f_cont[6]: xcovar=  nu.asfortranarray(xcovar)

    return avgloglikedata.contents.value

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
