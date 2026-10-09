TOBETESTED=si2
#
SHELL=/bin/sh
#
.PHONY: $(TOBETESTED)
#
all: $(TOBETESTED)
#
$(TOBETESTED):
	@echo "###########################################################################"
	@echo "running test $@ ..."
	cd $@; make all TEST=$@
#
