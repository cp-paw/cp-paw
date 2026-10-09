TOBETESTED=si2
#
SHELL=/bin/sh
#
.PHONY: $(TOBETESTED) skala
#
all: $(TOBETESTED)
#
$(TOBETESTED):
	@echo "###########################################################################"
	@echo "running test $@ ..."
	cd $@; make all TEST=$@

skala:
	@echo "###########################################################################"
	@echo "running optional Skala PAW integration test ..."
	cd skala_si2; make all
#
