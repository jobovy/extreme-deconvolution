SRC_FILES := $(filter-out src/main.c, $(wildcard src/*.c) $(wildcard *.c))
OBJECTS = $(addsuffix .o, $(notdir $(basename $(SRC_FILES))))
C = gcc
C_OPTS = -std=c99 -fPIC -shared -I. -Isrc
LINK_OPTS= -lm -lgsl -lgslcblas -lgomp

ExtremeDeconvolution.so: obj
	R CMD SHLIB -o ExtremeDeconvolution.so $(OBJECTS) $(LINK_OPTS)

obj: ; $(C) $(C_OPTS) -c $(SRC_FILES)

.PHONY: obj

clean:
	rm -f *.o *.so *.rds
