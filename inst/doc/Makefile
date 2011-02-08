# Makefile for LaTeX files

LATEX	= pdflatex
BIBTEX	= bibtex

SRC	= chemometrics-vignette.tex

all 	: chemometrics-vignette.pdf

%.pdf:	%.tex
	$(LATEX) $*
	$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*
	rm -f *.ps $*-* $*.dvi $*.log $*.bbl *.tex $*.aux $*.blg $*.log
	rm Rplots.pdf
