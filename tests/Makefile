TOBETESTED=unittests fulltests
#
SHELL=/bin/sh
#
.PHONY: $(TOBETESTED) profile
#
all: $(TOBETESTED)
#
$(TOBETESTED):
	@echo "###########################################################################"
	@echo "running $@ ..."
	cd $@; make all

profile:
	@echo "###########################################################################"
	@echo "running profiling cases ..."
	cd profile; make all
#
