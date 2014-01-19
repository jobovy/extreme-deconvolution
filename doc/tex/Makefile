all: ngerrors.pdf

%.pdf: %.ps
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true $<

%.ps: %.dvi
	dvips -t letter $< -o

%.dvi: %.tex
	latex $<
	latex $<
	latex $<
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && latex $< ) || echo noRerun "

#%.dvi: %.tex

.PHONY: clean spotless apj

clean:
	$(RM) *.log *.aux

spotless: clean
	$(RM) ngerrors.pdf