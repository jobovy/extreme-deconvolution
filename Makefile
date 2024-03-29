INSTALL_DIR=/usr/local/lib/
RM= /bin/rm -vf
R_EXEC=R
PYTHON=python
IDL=idl
ARCH=UNDEFINED
COVERAGE=0

ifeq ($(CC),)
  CC=gcc
endif

EDCFLAGS:= $(CFLAGS) -fopenmp -fcommon
EDLDFLAGS:= $(LDFLAGS) -fopenmp -fcommon

ifeq ($(ARCH),UNDEFINED)
	ARCH=$(shell uname -m)
endif

OS=$(shell uname -s)
ifeq ($(OS),Darwin)
	LIBEXT= dylib
	EDCFLAGS:= -arch $(ARCH) $(EDCFLAGS)
	EDLDFLAGS:= -arch $(ARCH) $(EDLDFLAGS)
	LINKOPTIONS:= -dynamiclib -Wl,-single_module
	RMDIR=rmdir
	ECHO=echo
	TRUE=TRUE
else
	LIBEXT= so
	LINKOPTIONS:= -shared
	RMDIR=rmdir -v
	ECHO=/bin/echo -e
	TRUE=true
endif
TARGETLIB= libextremedeconvolution.$(LIBEXT)

ifeq ($(COVERAGE),1)
	EDCFLAGS:= -O0 --coverage $(EDCFLAGS)
	EDLDFLAGS:= --coverage $(EDLDFLAGS)
endif

proj_gauss_mixtures_objects= src/bovy_randvec.o \
	src/calc_splitnmerge.o src/logsum.o src/minmax.o\
	src/normalize_row.o src/proj_EM.o src/proj_EM_step.o \
	src/proj_gauss_mixtures.o src/splitnmergegauss.o src/bovy_det.o

proj_gauss_main_objects= src/main.o src/parse_option.o src/read_data.o \
	src/read_IC.o src/read_till_sep.o src/write_model.o \
	src/cleanup.o

#
# The next targets are the main make targets: all, 
# extremedeconvolution (the executable), and 
# extremedeconvolution.so (the sharable object library)
#
all: build/$(TARGETLIB)

build:
	ls build 2>/dev/null ; \
	case "$$?" in \
	0);; \
	*) \
	mkdir build ;; \
	esac

# Build of main is currently broken, but this is *never* used, use the lib 
# instead
build/extremedeconvolution: $(proj_gauss_mixtures_objects) $(proj_gauss_main_objects) build
	$(CC) -o $@ \
	 $(EDCFLAGS) \
	 $(EDLDFLAGS) \
	 $(proj_gauss_mixtures_objects) \
	 $(proj_gauss_main_objects) \
	-lm -lgsl -lgslcblas -lgomp 2>/dev/null; \
	case "$$?" in \
	0);; \
	*) \
	$(CC) -o $@ \
	 $(EDCFLAGS) \
	 $(EDLDFLAGS) \
	 $(proj_gauss_mixtures_objects) \
	 $(proj_gauss_main_objects) \
	-lm -lgsl -lgslcblas;; \
	esac

build/$(TARGETLIB): $(proj_gauss_mixtures_objects) \
	src/proj_gauss_mixtures_IDL.o build
	$(CC) $(LINKOPTIONS) -o $@ \
	 $(EDLDFLAGS) \
	 $(proj_gauss_mixtures_objects) \
	 src/proj_gauss_mixtures_IDL.o \
	 -lm -lgsl -lgslcblas -lgomp 2>/dev/null; \
	case "$$?" in \
	0);; \
	*) \
	$(CC) $(LINKOPTIONS) -o $@ \
	 $(EDLDFLAGS) \
	 $(proj_gauss_mixtures_objects) \
	 src/proj_gauss_mixtures_IDL.o \
	 -lm -lgsl -lgslcblas ;; \
	esac

%.o: %.c
	$(CC) $(EDCFLAGS) -fpic -Wall -c $< -o $@ -I src/

#
# INSTALL THE IDL WRAPPER
#
install: build/$(TARGETLIB)
	cp $< $(INSTALL_DIR)$(TARGETLIB)

idlwrapper:
	(ls $(INSTALL_DIR)$(TARGETLIB) || ($(ECHO) "Cannot find library '$(INSTALL_DIR)$(TARGETLIB)', run again with 'INSTALL_DIR=' set to the directory you installed the library in" && exit -1))
	$(ECHO) 'result = CALL_EXTERNAL("$(INSTALL_DIR)$(TARGETLIB)", $$' > tmp
	cat pro/projected_gauss_mixtures_c.pro_1 tmp pro/projected_gauss_mixtures_c.pro_2 > pro/projected_gauss_mixtures_c.pro
	$(RM) tmp
	$(ECHO) 'Successfully installed IDL wrapper'

# INSTALL THE PYTHON WRAPPER
pywrapper:
	(ls $(INSTALL_DIR)$(TARGETLIB) || ($(ECHO) "Cannot find library '$(INSTALL_DIR)$(TARGETLIB)', run again with 'INSTALL_DIR=' set to the directory you installed the library in" && exit -1))
	(ls py/.extreme_deconvolution.py || mv py/extreme_deconvolution.py py/.extreme_deconvolution.py)
	sed "s#TEMPLATE_LIBRARY_PATH#'$(INSTALL_DIR)'#g" py/extreme_deconvolution_TEMPLATE.py > py/extreme_deconvolution.py
	((cd py && $(ECHO) 'import extreme_deconvolution' | $(PYTHON)) && $(ECHO) 'Successfully installed Python wrapper' || ($(ECHO) 'Something went wrong installing Python wrapper' && exit -1))

rpackage:
	cp src/*.h src/*.c r/src
	patch r/src/proj_gauss_mixtures_IDL.c < r/src/proj_gauss_mixtures_R.patch
	$(R_EXEC) CMD build r --no-manual
	$(R_EXEC) CMD check ExtremeDeconvolution_*.tar.gz --no-manual -o $(shell mktemp -d tmp.XXXX)
	rm -f r/src/*.h r/src/*.c

#
# TEST THE INSTALLATION
#
testidl:
	(cd examples && export IDL_PATH=$(IDL_PATH):../pro && $(ECHO) 'fit_TF' | $(IDL))
	(cd examples && ((diff TF.tex TF.out && $(ECHO) 'Ouput of test agrees with given solution') \
	|| $(ECHO) 'Output of test does not agree with given solution\nManually diff the TF.tex and TF.out (given solution) file'))

testpy:
	(cd py && $(PYTHON) extreme_deconvolution.py)


.PHONY: clean spotless

clean:
	$(RM) $(proj_gauss_mixtures_objects)
	$(RM) $(proj_gauss_main_objects)
	$(RM) src/proj_gauss_mixtures_IDL.o
	$(RM) $(subst src,r/src,$(proj_gauss_mixtures_objects))
	$(RM) $(subst src,r/src,$(proj_gauss_main_objects))
	$(RM) r/src/proj_gauss_mixtures_IDL.o
	$(RM) r/src/symbols.rds
	$(RM) ExtremeDeconvolution_*.tar.gz

spotless: clean rmbuild
	$(RM) src/*.~
	$(RM) pro/projected_gauss_mixtures.pro
	$(RM) py/*.pyc
	$(RM) pro/projected_gauss_mixtures_c.pro
	$(RM) examples/TF.ps examples/TF.tex
	$(RM) r/src/ExtremeDeconvolution.so
	(ls py/.extreme_deconvolution.py && mv py/.extreme_deconvolution.py py/extreme_deconvolution.py || $(TRUE))

rmbuild: build/
	$(RM) build/extremedeconvolution
	$(RM) build/$(TARGETLIB)
	($(RMDIR) build || $(ECHO) "Could not remove 'build/' directory, manually remove it")
