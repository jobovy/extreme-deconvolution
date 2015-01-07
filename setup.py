from distutils.core import setup, Extension
import os, os.path
import re

srcDir= 'src'
srcFiles= ['src/bovy_isfin.c','src/bovy_randvec.c','src/calc_splitnmerge.c',
		'src/logsum.c','src/minmax.c','src/normalize_row.c','src/proj_EM.c',
		'src/proj_EM_step.c','src/proj_gauss_mixtures.c',
		'src/splitnmergegauss.c','src/bovy_det.c',
		'src/proj_gauss_mixtures_IDL.c']


longDescription= "We present a general algorithm to infer a d-dimensional distribution function given a set of heterogeneous, noisy observations or samples. This algorithm reconstructs the error-deconvolved or 'underlying' distribution function common to all samples, even when the individual samples have unique error and missing-data properties. The underlying distribution is modeled as a mixture of Gaussians, which is completely general. Model parameters are chosen to optimize a justified, scalar objective function: the logarithm of the probability of the data under the error-convolved model, where the error convolution is different for each data point. Optimization is performed by an Expectation Maximization (EM) algorithm, extended by a regularization technique and 'split-and-merge' procedure. These extensions mitigate problems with singularities and local maxima, which are often encountered when using the EM algorithm to estimate Gaussian density mixtures."


setup(name='extreme-deconvolution',
      version='1.3',
      description='Density estimation using Gaussian mixtures in the presence of noisy, heterogeneous and incomplete data',
      author='Jo Bovy',
      author_email='jb2777@nyu.edu',
      license='BSD',
      long_description=longDescription,
      url='http://code.google.com/p/extreme-deconvolution/',
      package_dir = {'extreme_deconvolution': 'py'},
      ext_modules=[Extension('extreme_deconvolution._extreme_deconvolution', 
                             srcFiles,
                             libraries=['m','gsl','gslcblas','gomp'],
                             extra_compile_args=["-fopenmp"])],
      include_dirs=['src/'],
      packages=['extreme_deconvolution']
      )
