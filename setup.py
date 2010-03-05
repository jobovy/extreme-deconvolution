from distutils.core import setup, Extension
import os, os.path
import re

srcDir= 'src'
srcFiles= [os.path.join(srcDir,file) for file in os.listdir('src') if re.split(r'\.',file)[-1] == 'c']


longDescription= "We present a general algorithm to infer a d-dimensional distribution function given a set of heterogeneous, noisy observations or samples. This algorithm reconstructs the error-deconvolved or 'underlying' distribution function common to all samples, even when the individual samples have unique error and missing-data properties. The underlying distribution is modeled as a mixture of Gaussians, which is completely general. Model parameters are chosen to optimize a justified, scalar objective function: the logarithm of the probability of the data under the error-convolved model, where the error convolution is different for each data point. Optimization is performed by an Expectation Maximization (EM) algorithm, extended by a regularization technique and 'split-and-merge' procedure. These extensions mitigate problems with singularities and local maxima, which are often encountered when using the EM algorithm to estimate Gaussian density mixtures."


setup(name='extreme-deconvolution',
      version='1.2',
      description='Density estimation using Gaussian mixtures in the presence of noisy, heterogeneous and incomplete data',
      author='Jo Bovy',
      author_email='jb2777@nyu.edu',
      license='GPLv2',
      long_description=longDescription,
      url='http://code.google.com/p/extreme-deconvolution/',
      package_dir = {'': 'py'},
      py_modules=['extreme_deconvolution'],
      ext_modules=[Extension('extreme-deconvolution', srcFiles,libraries=['m','gsl','gslcblas'])],
      include_dirs=['src/'],
      data_files=[('doc',['extreme-deconvolution.pdf'])])
