from distutils.core import setup, Extension
import os, os.path
import re

srcDir= 'src'
srcFiles= [os.path.join(srcDir,file) for file in os.listdir('src') if re.split(r'\.',file)[-1] == 'c']

setup(name='extreme-deconvolution',
      version='1.1',
      description='Density estimation using Gaussian mixtures in the presence of noisy, heterogeneous and incomplete data',
      author='Jo Bovy',
      author_email='jb2777@nyu.edu',
      url='http://code.google.com/p/extreme-deconvolution/',
      package_dir = {'': 'py'},
      py_modules=['extreme_deconvolution'],
      ext_modules=[Extension('extreme-deconvolution', srcFiles)],
      include_dirs=['src/'])
