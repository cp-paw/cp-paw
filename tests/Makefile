TOBETESTED=unittests fulltests
#
SHELL=/bin/sh
#
.PHONY: $(TOBETESTED)
#
all: $(TOBETESTED)
#
$(TOBETESTED):
	@echo "###########################################################################"
	@echo "running $@ ..."
	cd $@; make all
#
