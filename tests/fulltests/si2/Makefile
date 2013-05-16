SHELL=/bin/sh
#
###########################################################################
## replaces paw.f90 with paw_mod_for_unittests.f90 (no MAIN PRORGAM)     ##
## compiles ${TESTS}.f90 and runs paw_${TEST}.x                          ##
###########################################################################
all:
	@echo "running paw..."
	${PAWX} ${TEST}.cntl 2> err.log > out.log
	@echo "analyse results ..."
	@${SHELL} analyse.sh 2> err.log > out.log; if [ $$? -eq 0 ] ; then echo "test ${TEST} passed."; else echo "test ${TEST} failed. See `pwd`/out.log and `pwd`/err.log."; exit 1 ; fi
	
