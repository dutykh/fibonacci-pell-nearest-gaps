# Rational reflections and palindromic words in the Cohn monoid
#
# Authors:
#   Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#   and Technology, Abu Dhabi, UAE)
#   Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery,
#   France)
#
MAIN :=DD-LV-Palindromic-Reflection-Rigidity
PYTHON := python3
LATEXMK := latexmk
SOURCES := $(MAIN).tex references.bib $(wildcard sections/*.tex) $(wildcard figures/*.tex)

.PHONY: all help rebuild check checks release clean
all: $(MAIN).pdf

$(MAIN).pdf: $(SOURCES)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PYTHON) -B tools/check_build.py $(MAIN)

help:
	@echo 'Targets:'
	@echo '  all      build $(MAIN).pdf (default)'
	@echo '  rebuild  force a full rebuild with the same strict checks'
	@echo '  check    alias for rebuild; fails on any LaTeX or BibTeX warning'
	@echo '  checks   run the exact computations of supplement/'
	@echo '  release  rebuild the manuscript, run the computations, then clean'
	@echo '  clean    remove LaTeX intermediates, keep the PDF'

rebuild:
	$(LATEXMK) -C $(MAIN).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PYTHON) -B tools/check_build.py $(MAIN)

check: rebuild

checks:
	$(PYTHON) -B supplement/run_all.py

release: rebuild checks
	$(MAKE) clean

clean:
	$(LATEXMK) -c $(MAIN).tex
	$(RM) $(MAIN).bbl
