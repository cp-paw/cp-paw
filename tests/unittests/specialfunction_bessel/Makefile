SHELL=/bin/sh
#
###########################################################################
## replaces paw.f90 with paw_mod_for_unittests.f90 (no MAIN PRORGAM)     ##
## compiles ${TESTS}.f90 and runs paw_${TEST}.x                          ##
###########################################################################
all:
	@#preprocess ${UNITTESTS}.f90
	@${F90PP} ${TEST}.${FEXT} >${OBJDIR}/paw_unittest_${TEST}_d.${FEXT}
	@#compile ${UNITTESTS}.f90
	@${COMPILE} -o ${OBJDIR}/paw_unittest_${TEST}_d.o ${OBJDIR}/paw_unittest_${TEST}_d.${FEXT}
	@#link binary
	@${LINK} -o paw_unittest_${TEST}.x $(subst paw.o,paw_unittest_${TEST}_d.o,${OBJECTS}) $(subst paw.o,paw_unittest_${TEST}_d.o,${LOBJS}) ${OBJVERSION} ${LOBJSF77} ${OBJDIR}/version_writeparmfile.o ${MPE} ${LIBS}
	@#run unittest
	@./paw_unittest_${TEST}.x
