RM= /bin/rm -vf

proj_gauss_mixtures_objects= src/bovy_isfin.o src/bovy_randvec.o \
	src/calc_splitnmerge.o src/logsum.o src/minmax.o\
	src/normalize_row.o src/proj_EM.o src/proj_EM_step.o \
	src/proj_gauss_mixtures.o src/splitnmergegauss.o src/bovy_det.o
#	src/calc_qstarij.o

proj_gauss_main_objects= src/main.o src/parse_option.o src/read_data.o \
		src/read_IC.o src/read_till_sep.o src/write_model.o \
		src/cleanup.o

#
# The next targets are the main make targets: all, 
# proj_gauss_mixtures (the executable), and 
# proj_gauss_mixtures.so (the sharable object library)
#
all: proj_gauss_mixtures proj_gauss_mixtures.so

proj_gauss_mixtures: $(proj_gauss_mixtures_objects) $(proj_gauss_main_objects)
	gcc -o proj_gauss_mixtures -lm -lgsl -lgslcblas \
	$(proj_gauss_mixtures_objects) $(proj_gauss_main_objects)

proj_gauss_mixtures.so: $(proj_gauss_mixtures_objects) \
			src/proj_gauss_mixtures_IDL.o
	gcc -shared -o proj_gauss_mixtures.so -lm -lgsl -lgslcblas \
	$(proj_gauss_mixtures_objects) src/proj_gauss_mixtures_IDL.o

#
# proj_gauss_main_objects
#
src/main.o: src/main.c
	gcc -fpic -Wall -c src/main.c -o src/main.o -I src/

src/read_till_sep.o: src/read_till_sep.c
	gcc -fpic -Wall -c src/read_till_sep.c -o src/read_till_sep.o -I src/

src/read_data.o: src/read_data.c
	gcc -fpic -Wall -c src/read_data.c -o src/read_data.o -I src/

src/read_IC.o: src/read_IC.c
	gcc -fpic -Wall -c src/read_IC.c -o src/read_IC.o -I src/

src/parse_option.o: src/parse_option.c
	gcc -fpic -Wall -c src/parse_option.c -o src/parse_option.o -I src/

src/write_model.o: src/write_model.c
	gcc -fpic -Wall -c src/write_model.c -o src/write_model.o -I src/

src/cleanup.o: src/cleanup.c
	gcc -fpic -Wall -c src/cleanup.c -o src/cleanup.o -I src/


#
# proj_gauss_mixtures_objects
#
src/bovy_isfin.o: src/bovy_isfin.c
	gcc -fpic -Wall -c src/bovy_isfin.c -o src/bovy_isfin.o

src/bovy_randvec.o: src/bovy_randvec.c
	gcc -fpic -Wall -c src/bovy_randvec.c -o src/bovy_randvec.o -I src/

src/bovy_det.o: src/bovy_det.c
	gcc -fpic -Wall -c src/bovy_det.c -o src/bovy_det.o -I src/

#src/calc_qstarij.o: src/calc_qstarij.c
#	gcc -fpic -Wall -c src/calc_qstarij.c -o src/calc_qstarij.o \
#	-I src/

src/calc_splitnmerge.o: src/calc_splitnmerge.c
	gcc -fpic -Wall -c src/calc_splitnmerge.c -o src/calc_splitnmerge.o \
	-I src/

src/splitnmergegauss.o: src/splitnmergegauss.c
	gcc -fpic -Wall -c src/splitnmergegauss.c -o src/splitnmergegauss.o \
	-I src/

src/logsum.o: src/logsum.c
	gcc -fpic -Wall -c src/logsum.c -o src/logsum.o -I src/

src/minmax.o: src/minmax.c
	gcc -fpic -Wall -c src/minmax.c -o src/minmax.o -I src/

src/normalize_row.o: src/normalize_row.c
	gcc -fpic -Wall -c src/normalize_row.c -o src/normalize_row.o -I src/

src/proj_EM.o: src/proj_EM.c
	gcc -fpic -Wall -c src/proj_EM.c -o src/proj_EM.o -I src/

src/proj_EM_step.o: src/proj_EM_step.c
	gcc -fpic -Wall -c src/proj_EM_step.c -o src/proj_EM_step.o -I src/

src/proj_gauss_mixtures.o: src/proj_gauss_mixtures.c
	gcc -fpic -Wall -c src/proj_gauss_mixtures.c -o \
	src/proj_gauss_mixtures.o -I src/


#
# IDL wrapper routine
#
src/proj_gauss_mixtures_IDL.o: src/proj_gauss_mixtures_IDL.c
	gcc -fpic -Wall -c src/proj_gauss_mixtures_IDL.c -o \
	src/proj_gauss_mixtures_IDL.o -I src/


#
# INSTALL THE IDL WRAPPER
#
install:
	echo 'result = CALL_EXTERNAL("$(PWD)/proj_gauss_mixtures.so", $$' > tmp
	cat pro/projected_gauss_mixtures_c.pro_1 tmp pro/projected_gauss_mixtures_c.pro_2 > pro/projected_gauss_mixtures_c.pro
	$(RM) tmp


#
# TEST THE INSTALLATION
#
test:
	(cd examples && echo 'fit_TF' | idl)
	(cd examples && ((diff TF.tex TF.out && echo 'Ouput of test agrees with given solution') \
	|| echo -e 'Output of test does not agree with given solution\nManually diff the TF.tex and TF.out (given solution) file'))


.PHONY: clean spotless

clean:
	$(RM) $(proj_gauss_mixtures_objects)
	$(RM) $(proj_gauss_main_objects)
	$(RM) src/proj_gauss_mixtures_IDL.o

spotless: clean
	$(RM) proj_gauss_mixtures.so
	$(RM) proj_gauss_mixtures
	$(RM) src/*.~
