SRC_FILES := $(filter-out main.c, $(wildcard *.c))
OBJECTS    = $(addsuffix .o, $(notdir $(basename $(SRC_FILES))))
PKG_CFLAGS = -std=c99 -fPIC -fopenmp -shared -I. -Isrc
PKG_CLIBS  = -lm -lgsl -lgslcblas -lgomp

ExtremeDeconvolution.so: obj
	R CMD SHLIB -o ExtremeDeconvolution.so $(OBJECTS) $(PKG_CLIBS)

obj: ; $(CC) $(CFLAGS) $(PKG_CFLAGS) -c $(SRC_FILES)

.PHONY: obj

clean:
	rm -f *.o *.so *.rds
