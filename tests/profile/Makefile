TOBETESTED=si64
#
SHELL=/bin/sh
#
.PHONY: all clean $(TOBETESTED)
#
all: $(TOBETESTED)
#
$(TOBETESTED):
	@echo "###########################################################################"
	@echo "running profile case $@ ..."
	cd $@; make all TEST=$@

clean:
	@for dir in $(TOBETESTED); do \
	  cd $$dir; make clean; cd ..; \
	done
