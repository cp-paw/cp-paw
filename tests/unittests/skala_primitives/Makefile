SHELL=/bin/sh

all:
	@${F90PP} <${TEST}.${FEXT} >${OBJDIR}/paw_unittest_${TEST}_d.${FEXT}
	@${COMPILE} -o ${OBJDIR}/paw_unittest_${TEST}_d.o ${OBJDIR}/paw_unittest_${TEST}_d.${FEXT}
	@${LINK} -o paw_unittest_${TEST}.x $(subst paw.o,paw_unittest_${TEST}_d.o,${OBJECTS}) $(subst paw.o,paw_unittest_${TEST}_d.o,${LOBJS}) ${OBJVERSION} ${LOBJSF77} ${OBJDIR}/version_writeparmfile.o ${MPE} ${LIBS}
	@./paw_unittest_${TEST}.x
