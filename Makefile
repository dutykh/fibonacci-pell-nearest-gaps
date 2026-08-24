# Fibonacci-Pell nearest gaps: an all-exponent classification and
# quadratic-unit orbit rigidity
#
# Authors:
#   Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#   and Technology, Abu Dhabi, UAE)
#   Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery,
#   France)
#
LATEXMK := latexmk
PYTHON  := python3
MAIN := DD-LV-Fibonacci-Pell-Gaps
SOURCES := $(MAIN).tex references.bib $(wildcard sections/*.tex)

.PHONY: all help rebuild clean distclean check certificates release

# The build is gated on a warning-free log. Two warnings are whitelisted by
# exact text, and only those two:
#   * the pdfTeX font-expansion notice, which is emitted before first use and
#     has no effect on the output;
#   * amsplain's missing-pages warning for PomeoBravo2024, which is online-first
#     with no volume, issue or pages assigned (verified against its Crossref
#     record). Pagination is not invented to silence it.
define compile_and_check
$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
@test -s $(MAIN).log
@test -s $(MAIN).blg
@! grep -E "LaTeX (Font )?Warning:|Package .* Warning:|Class .* Warning:|Underfull \\\\[hv\\\\]box|Overfull \\\\[hv\\\\]box|Missing character:" $(MAIN).log
@! grep -E "pdfTeX warning" $(MAIN).log | grep -Fv "pdfTeX warning (font expansion): font should be expanded before its first use"
@! grep -E "multiply defined" $(MAIN).log
@! grep -E "^Warning--" $(MAIN).blg | grep -Fv "missing pages in PomeoBravo2024"
@! grep -n -P '\\(?:leq?|geq?)(?!slant|[A-Za-z])' $(MAIN).tex sections/*.tex
$(LATEXMK) -c $(MAIN).tex
$(RM) $(MAIN).bbl
endef

all: $(MAIN).pdf

help:
	@echo 'Targets:'
	@echo '  all           build $(MAIN).pdf (default)'
	@echo '  rebuild       force a full rebuild with the same strict checks'
	@echo '  check         alias for rebuild; fails on any LaTeX or BibTeX warning'
	@echo '  certificates  run every exact Python certificate in supplement/'
	@echo '  release       rebuild the manuscript, then run every certificate'
	@echo '  clean         remove LaTeX intermediates, keep the PDF'
	@echo '  distclean     remove intermediates AND the tracked PDF'

$(MAIN).pdf: $(SOURCES)
	$(compile_and_check)

rebuild:
	$(compile_and_check)

check: rebuild

certificates:
	$(PYTHON) -B supplement/run_all.py

release: rebuild certificates

clean:
	$(LATEXMK) -c $(MAIN).tex
	$(RM) $(MAIN).bbl

distclean:
	$(LATEXMK) -C $(MAIN).tex
