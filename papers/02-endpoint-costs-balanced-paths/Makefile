# Endpoint costs and balanced near-minimizers under Gaussian reflections
#
# Authors:
#   Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#   and Technology, Abu Dhabi, UAE)
#   Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery,
#   France)
#
MAIN := DD-LV-Endpoint-Costs-Balanced-Paths
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
	@echo '  checks   run every exact Python program in supplement/'
	@echo '  release  rebuild the manuscript, run every program, then clean'
	@echo '  clean    remove LaTeX intermediates, keep the PDF'

rebuild:
	$(LATEXMK) -C $(MAIN).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PYTHON) -B tools/check_build.py $(MAIN)

check: rebuild

checks:
	$(PYTHON) -B supplement/verify.py
	$(PYTHON) -B supplement/verify_independent.py
	$(PYTHON) -B supplement/verify_dictionary.py

release: rebuild checks
	$(MAKE) clean

clean:
	$(LATEXMK) -c $(MAIN).tex
	$(RM) $(MAIN).bbl
