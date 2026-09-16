# Around the Markoff uniqueness conjecture: a series of manuscripts and their
# reproduction material.
#
# Authors:
#   Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
#   and Technology, Abu Dhabi, UAE)
#   Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery,
#   France)
#
# This file builds nothing itself. It discovers the manuscript directories
# below papers/ and forwards each target to the Makefile each one carries, so
# a new manuscript joins the series without an edit here. Every manuscript
# directory is expected to provide the five common targets
#
#   all      build the manuscript PDF
#   check    force a full rebuild under that manuscript's strict gates
#   checks   run the reproduction programs of its supplement
#   release  rebuild and then run the reproduction programs
#   clean    remove the build intermediates
#
# and may provide further targets of its own, reachable with
# make -C papers/<directory> <target>.

PAPERS := $(sort $(notdir $(patsubst %/,%,$(dir $(wildcard papers/*/Makefile)))))

.PHONY: all help list check checks release clean $(PAPERS)

all: $(addprefix build-,$(PAPERS))

help:
	@echo 'Series targets, applied to every manuscript below papers/:'
	@echo '  all        build every manuscript PDF (default)'
	@echo '  check      force a full rebuild of every manuscript'
	@echo '  checks     run the reproduction programs of every supplement'
	@echo '  release    rebuild every manuscript, then run every supplement'
	@echo '  clean      remove build intermediates everywhere'
	@echo '  list       list the manuscript directories that were discovered'
	@echo ''
	@echo 'One manuscript at a time:'
	@echo '  make <directory>                build that manuscript'
	@echo '  make -C papers/<directory> help its own targets'
	@echo ''
	@echo 'Discovered manuscripts:'
	@$(foreach p,$(PAPERS),echo '  $(p)';)

list:
	@$(foreach p,$(PAPERS),echo '$(p)';)

check: $(addprefix check-,$(PAPERS))
checks: $(addprefix checks-,$(PAPERS))
release: $(addprefix release-,$(PAPERS))
clean: $(addprefix clean-,$(PAPERS))

# For each discovered manuscript, a bare target that builds it and one
# forwarding target per verb of the common contract.
define paper_rules
.PHONY: build-$(1) check-$(1) checks-$(1) release-$(1) clean-$(1)
$(1): build-$(1)
build-$(1):
	$$(MAKE) -C papers/$(1) all
check-$(1):
	$$(MAKE) -C papers/$(1) check
checks-$(1):
	$$(MAKE) -C papers/$(1) checks
release-$(1):
	$$(MAKE) -C papers/$(1) release
clean-$(1):
	$$(MAKE) -C papers/$(1) clean
endef

$(foreach p,$(PAPERS),$(eval $(call paper_rules,$(p))))
