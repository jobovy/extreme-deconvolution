Extreme Deconvolution (XD)
===========================

**Density estimation using Gaussian mixtures in the presence of noisy, heterogeneous and incomplete data**

.. image:: https://github.com/jobovy/extreme-deconvolution/actions/workflows/build.yml/badge.svg
   :target: https://github.com/jobovy/extreme-deconvolution/actions/workflows/build.yml
   
Summary
---------

Extreme-deconvolution (XD) is a general algorithm to infer a d-dimensional distribution function from a set of heterogeneous, noisy observations or samples. It is fast, flexible, and treats the data's individual uncertainties properly, to get the best description possible of the underlying distribution. It performs well over the full range of density estimation, from small data sets with only tens of samples per dimension, to large data sets with millions of data points.

The extreme-deconvolution algorithm is available here as a dynamic C-library that your programs can link to, or through Python, R, or IDL wrappers that allow you to call the fast underlying C-code in your high-level applications with minimal overhead.

News
------

* **2016/11/04**: Check out the `XDGMM <https://github.com/tholoien/XDGMM>`__ wrapper class, which allows fitting, model selection, resampling, and conditioning of the Gaussian mixture model. Also check out the `astroML <http://www.astroml.org/index.html>`__ alternative implementation of the XD algorithm (Python only).

* **2015/08/17**: An interface in R was added by `@gaow <https://github.com/gaow>`__. See the updated "Installation notes" and "Usage notes" below.

* **2014/01/17**: Code was migrated to github; previous news refers to the googlecode version.

* **2011/06/12**: Version 1.3 released. This version incorporates OpenMP multiprocessing support and other improvements since version 1.2 (see other _News_ below).  

* **2011/01/29**: The code now has multiprocessing support using OpenMP and is automatically compiled with OpenMP. Use the environment variable OMP_NUM_THREADS to control the number of threads used. Only available in the `trunk <http://code.google.com/p/extreme-deconvolution/source/browse/trunk>`__ at this time. If you want to compile the code without OpenMP, remove the -fopenmp references in the Makefile and copy the src/_omp.h file to src/omp.h.

* **2010/06/11**: Some `addons <https://github.com/jobovy/extreme-deconvolution/tree/master/addons>`__ were added: Use these small IDL programs to use the output from the core extreme-deconvolution code (e.g., calculate uncertainties on the best-fit parameters, calculate membership probabilities, perform a KS test with these membership probabilities). Download them `here <http://extreme-deconvolution.googlecode.com/files/extreme-deconvolution-addons_1.0.tar.gz>`__ or check-out the code.

* **2010/04/01**: The "weight" option was added to give the data points different weights in the calculation of the log likelihood (log likelihood = weight `*` log p(data|model) ). Only available in the `trunk <http://code.google.com/p/extreme-deconvolution/source/browse/trunk>`__ at this time.

* **2010/03/05**: Version 1.2 released. The "projection" matrices are now optional and the data-uncertainties can be diagonal (i.e., uncorrelated).

Requirements
------------

`GSL <http://www.gnu.org/software/gsl/>`__: The GNU Scientific Library. See for example `this page <https://galpy.readthedocs.io/en/latest/installation.html#how-do-i-install-the-gsl>`__ for information on how to install the GSL on your system.

Get the latest version
----------------------

Get the latest version by checking out the git repository:

``git clone https://github.com/jobovy/extreme-deconvolution.git``

or download it using the big green button above the file listing above. 

(note that downloading and installing the latest released version under the 'releases' tab is *not* recommended, as this version is out-of-date; please install the latest master version instead)

Installation
------------

To compile the code, navigate to the directory where you downloaded the code and do::

   make

To install the library do::

   sudo make install

or::

   make install INSTALL_DIR=/path/to/install/dir/


To install the IDL wrapper do::

   make idlwrapper

Add INSTALL_DIR=/path/to/install/dir/ if you used this to install the library


To install the Python wrapper do::

   make pywrapper

Add INSTALL_DIR=/path/to/install/dir/ if you used this to install the library. Remember to add the ``py/`` directory to your `PYTHONPATH` to use the code in Python.

To install the R package do::

   make rpackage
   R CMD INSTALL ExtremeDeconvolution_1.3.tar.gz

Fix the version number as needed. Note that options for compiling
packages in R are specified through the ``Makevars`` file, which
should typically be located at ``~/.R/Makevars``. For example, if you
need to override the default C compiler to ``gcc-4.9``, you would add
line ``CC=gcc-4.9`` to the Makevars file before building the
package. (You also need to make sure that the proper ``CC`` is set in
the main ``Makefile`` as well.)  For more details on customzing R package
installation, see `here
<http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Customizing-package-compilation>`__.
Alternatively, you may find that it is more convenient to use the
`install.packages()` function in R to install the package. In that
case, replace the second step (``R CMD INSTALL ...``) with the
following call within your R environment ::

   install.packages(pkgs = "ExtremeDeconvolution_1.3.tar.gz",repos = NULL)

This assumes that the R working directory is the same as the root of
this git repository.
   
To test whether the code and the python wrapper is working do::

   make testpy

To test whether the code and the IDL wrapper is working do (requires IDL and the IDL-wrapper to be installed)::

   make testidl

Clean up intermediate files::

   make clean

Usage
------

Examples of use of the code are in the IDL example code in `<examples/fit_tf.pro>`__ and in the python doctest in `<py/extreme_deconvolution.py>`__.

In python you would typically do something like::

   from extreme_deconvolution import extreme_deconvolution
   #Set up your arrays: ydata has the data, ycovar the uncertainty covariances
   #initamp, initmean, and initcovar are initial guesses
   #get help on their shapes and other options using
   ?extreme_deconvolution
   #Run the code
   extreme_deconvolution(ydata,ycovar,initamp,initmean,initcovar)
   #initamp, initmean, and initcovar are now updated to their best fit values


In IDL this becomes::

   ;;Set up arrays and the number of Gaussians
   K=2 ;;K Gaussians
   ;;Run the code
   projected_gauss_mixtures_c, K, ydata, ycovar, initamp, initmean, initcovar, /quiet
   ;;initamp, initmean, and initcovar are now updated to their best fit values


In R::

   library(ExtremeDeconvolution)
   ?extreme_deconvolution


Installation FAQ
-----------------

* *`make` returns "file was built for unsupported file format which is not the architecture being linked (i386)" errors (or x86_64)*

  XD is trying to compile as a 32 (or 64) bit library while your GSL or OpenMP libraries were compiled as 64 (or 32) bit libraries. You can force XD to compile as a particular architecture by adding the ARCH option to make, e.g.::

     make ARCH=x86_64


* *I do not have/want OpenMP*

  You can disable OpenMP support by removing the -fopenmp and -lgomp references in the Makefile. If you are installing the R package, also remove the -fopenmp and -lgomp references in r/src/Makefile.

* *Problems with clang*

  On Macs with OS X >= 10.9, gcc is no longer the default compiler, which is instead clang (although confusingly, gcc points to clang!). Clang does not have support for OpenMP (yet) and the code will therefore only run on a single CPU. To use the OpenMP parallelized version of the code, install gcc yourself and make sure that the Makefile is using it (using the CC variable). One recommended option on a Mac is to install gcc with openmp using `Homebrew <http://brew.sh>`__; e.g., after installing Homebrew on your Mac, run ``brew install gcc@4.9``, then set ``CC=gcc-4.9`` in the Makefile for this repository.

Acknowledgments
-----------------

Thanks to Gao Wang and Peter Carbonetto for the R interface and
Daniela Carollo, Joe Hennawi, Sergey Koposov, and Leonidas Moustakas
for bug reports and fixes.

Acknowledging extreme-deconvolution
------------------------------------

The algorithm that the code implements is described in the paper *Extreme deconvolution: inferring complete distribution functions from noisy, heterogeneous and incomplete observations*; a copy of the latest draft of this paper is included in the "doc/" directory of the repository or source archive. If you use the code, please cite this paper, e.g.::

    Extreme deconvolution: inferring complete distribution functions from noisy, heterogeneous and incomplete observations
    Jo Bovy, David W. Hogg, & Sam T. Roweis, Ann. Appl. Stat. 5, 2B, 1657 (2011)


Extreme-deconvolution in action
--------------------------------

A good place to find examples is the `citations to the extreme-deconvolution paper <http://adsabs.harvard.edu/abs/2011AnApS...5.1657B>`__. The code is also used in a variety of fields outside of astronomy.
