# Euclidean minimality of Markov representations on thirteen Farey rows
#
# Authors:
#   Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#   and Technology, Abu Dhabi, UAE)
#   Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery,
#   France)
#
# The strict build gate of this manuscript is part of its supplement:
# supplement/code/check_build.py.
MAIN := DD_LV_Thirteen_Farey_Rows
PYTHON := python3
LATEXMK := latexmk
SOURCES := $(MAIN).tex references.bib $(wildcard sections/*.tex) $(wildcard figures/*.tex)

.PHONY: all help rebuild check checks figures release clean
all: $(MAIN).pdf

$(MAIN).pdf: $(SOURCES)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PYTHON) -B supplement/code/check_build.py $(MAIN)

help:
	@echo 'Targets:'
	@echo '  all      build $(MAIN).pdf (default)'
	@echo '  rebuild  force a full rebuild with the same strict checks'
	@echo '  check    alias for rebuild; fails on any LaTeX or BibTeX warning'
	@echo '  checks   run the thirteen exact checks of supplement/'
	@echo '  figures  rewrite the TikZ code of the figures from exact data'
	@echo '  release  rebuild the manuscript, run the checks, then clean'
	@echo '  clean    remove LaTeX intermediates, keep the PDF'

rebuild:
	$(LATEXMK) -C $(MAIN).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PYTHON) -B supplement/code/check_build.py $(MAIN)

check: rebuild

checks:
	$(PYTHON) -B supplement/run_all.py

figures:
	$(PYTHON) -B supplement/code/make_figures.py

release: rebuild checks
	$(MAKE) clean

clean:
	$(LATEXMK) -c $(MAIN).tex
	$(RM) $(MAIN).bbl
