#########################################################################
#########################################################################
##                                                                     ##
##  make an paw executable and a paw library using gmake               ##
##  ====================================================               ##
##                                                                     ##
##                                                                     ##
##  The following environment variables must be set before             ##
##  invoking paw.mk                                                    ##
##  SRCDIR    holds the *.f90 files                                    ##
##  OBJDIR    holds the *.o files                                      ##
##  PAWX      is the name of the executable                            ##
##  F90PP     fortran90 preprocessor (self-made)                       ##
##  COMPILE   compile command including all options                    ##
##  COMPILE77 compile command for fortran77 including all options      ##
##  FEXT      extension required by the compiler for source files      ##
##  LINK      link command including all options                       ##
##  PARALLEL  PARALLEL=T compiles the parallel version of the program  ##
##  ALLDEP    ALLDEP=T uses all module dependencies                    ##
##  OPTRM     comand executed on intermediate source files (remove or touch)##
#########################################################################
#
##########################################################################
##  define parameters                                                   ##
##########################################################################
#
ifeq ($(PARALLEL),T)   #select parallel versus scalar paw_mpelib.f90
  MODE = _p
  MPE = ${OBJDIR}paw_mpelib$(MODE).o
else
  MODE = _s
  MPE = ${OBJDIR}paw_mpelib$(MODE).o
endif
#          leave the hatch right after ...a in the definition of PAWLIB
PAWLIB   = ${OBJDIR}libpaw.a#         paw library
#
##########################################################################
##  define object files kept in the library                             ##
##########################################################################
#        $(PAWLIB)(paw_mpelib$(MODE).o) \
#
LOBJS = \
        ${OBJDIR}paw_trace.o \
        ${OBJDIR}paw_error.o \
        ${OBJDIR}paw_filehandler.o \
        ${OBJDIR}paw_clock.o \
        ${OBJDIR}paw_lock.o \
        ${OBJDIR}paw_timing.o \
        ${OBJDIR}paw_linkedlist.o \
        ${OBJDIR}paw_strings.o \
        ${OBJDIR}paw_dft.o \
        ${OBJDIR}paw_dftaddendum.o \
        ${OBJDIR}paw_constants.o \
        ${OBJDIR}paw_spherical.o \
        ${OBJDIR}paw_generalpurpose.o \
        ${OBJDIR}paw_report.o \
        ${OBJDIR}paw_periodictable.o \
        ${OBJDIR}paw_radial.o \
        ${OBJDIR}paw_usage.o \
        ${OBJDIR}paw_selftest.o \
        ${OBJDIR}paw_strcio.o \
        ${OBJDIR}paw_cell.o \
        ${OBJDIR}paw_formats.o \
        ${OBJDIR}paw_library.o \
        ${OBJDIR}paw_polynom.o 
#
##########################################################################
##  define other objects                                                ##
##########################################################################
#
OBJECTS1 = \
           ${OBJDIR}paw_thermostat.o \
           ${OBJDIR}paw_isolate.o \
           ${OBJDIR}paw_assist.o \
           ${OBJDIR}paw_ldaplusu.o \
           ${OBJDIR}paw_lists.o \
           ${OBJDIR}paw_constraints.o \
           ${OBJDIR}paw_fft.o \
#
OBJECTS2 = \
           ${OBJDIR}paw_augmentation.o \
           ${OBJDIR}paw_atoms.o \
           ${OBJDIR}paw.o \
           ${OBJDIR}paw_efg.o \
           ${OBJDIR}paw_ioroutines.o \
           ${OBJDIR}paw_iotra.o \
           ${OBJDIR}paw_ionew.o \
           ${OBJDIR}paw_qmmm.o \
           ${OBJDIR}paw_classical.o \
           ${OBJDIR}paw_continuum.o \
           ${OBJDIR}paw_warmup.o \
#  the following is for the QMMM of Woo and Margl
#          ${OBJDIR}mm_paw_modules.o \
#          ${OBJDIR}mm_paw_core_mm.o \
#          ${OBJDIR}mm_paw_interface.o
#          ${OBJDIR}paw_md
#
OBJECTS3 = \
           ${OBJDIR}paw_setups.o \
           ${OBJDIR}paw_potential.o \
           ${OBJDIR}paw_occupations.o \
           ${OBJDIR}paw_pairpotential.o \
           ${OBJDIR}paw_graphics.o \
           ${OBJDIR}paw_waves.o \
           ${OBJDIR}paw_pdos.o \
           ${OBJDIR}paw_kpoints.o \
           ${OBJDIR}paw_vext.o
#          ${OBJDIR}paw_optic.o \
#
OBJECTS= ${OBJECTS1} ${OBJECTS2} ${OBJECTS3} 
# 
#
#########################################################################
##  create preprocessor                                                ##
#########################################################################
F90PP_TMPLTS=${PAWDIR}F90PP/f90pp_tmplts
${F90PP_TMPLTS}.x : ${F90PP_TMPLTS}.f
	${F90PP} -notemplates $< >${F90PP_TMPLTS}_d.${FEXT}
	${COMPILE} -o ${F90PP_TMPLTS}.o ${F90PP_TMPLTS}_d.${FEXT}
	${LINK} -o ${F90PP_TMPLTS}.x ${F90PP_TMPLTS}.o
#
#########################################################################
##  link files and produce executable PAWX                             ##
#########################################################################
#
${PAWX} : ${LOBJS} ${OBJECTS} ${MPE} ${OBJDIR}liblapack.a ${OBJDIR}libblas.a
	${LINK} -o ${PAWX} ${OBJECTS} ${LOBJS} ${MPE} ${LLIBS}
	ar -ru $(PAWLIB) $(LOBJS) ${OBJDIR}paw_mpelib_s.o
#
#${OBJDIR}mm_paw_f77.o 
#########################################################################
##  compile code files independent of parameter                        ##
##  and place the objects in OBJDIR                                    ##
#########################################################################
#
#__special rule for parallel interface because of preprocessor___________
${MPE} : ${OBJDIR}%$(MODE).o: ${SRCDIR}%.f
	${F90PP} $< >${OBJDIR}$*$(MODE)_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}$*$(MODE)_d.${FEXT}
	${OPTRM} ${OBJDIR}$*$(MODE)_d.${FEXT}
#
#__compile lapack routines________________________________________________
${OBJDIR}liblapack.a : ${SRCDIR}paw_lapack.f
	${COMPILE77} -o ${OBJDIR}paw_lapack.o ${SRCDIR}paw_lapack.f
	ar -ru ${OBJDIR}liblapack.a ${OBJDIR}paw_lapack.o
	${OPTRM} ${OBJDIR}paw_lapack.o
#
#__compile blas routines________________________________________________
${OBJDIR}libblas.a : ${SRCDIR}paw_blas.f
	${COMPILE77} -o ${OBJDIR}paw_blas.o ${SRCDIR}paw_blas.f
	ar -ru ${OBJDIR}libblas.a ${OBJDIR}paw_blas.o
	${OPTRM} ${OBJDIR}paw_blas.o
#
#__general rule for object files_________________________________________
$(OBJECTS) $(LOBJS): ${OBJDIR}%.o: ${SRCDIR}%.f
	${F90PP} $< >${OBJDIR}$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}$*_d.${FEXT}
	${OPTRM} ${OBJDIR}$*_d.${FEXT}
#
#########################################################################
##  include dependencies through mod files                             ##
#########################################################################
ifeq ($(ALLDEP),T)
 ${OBJDIR}paw_filehandler.o:   ${OBJDIR}paw_strings.o
 ${OBJDIR}paw_linkedlist.o:    ${OBJDIR}paw_strings.o
 ${OBJDIR}paw_periodictable.o: ${OBJDIR}paw_strings.o
 ${OBJDIR}paw_ioroutines.o:    ${OBJDIR}paw_strings.o
#
 ${OBJDIR}paw_timing.o:        ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_linkedlist.o:    ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_report.o:        ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_waves.o:         ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_augmentation.o:  ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw.o:               ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_constraints.o:   ${OBJDIR}paw_mpelib${MODE}.o 
 ${OBJDIR}paw_potential.o:     ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_pairpotential.o: ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_isolate.o:       ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_qmmm.o:          ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_thermostat.o:    ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_occupations.o:   ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_ionew.o:         ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_atoms.o:         ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_classical.o:     ${OBJDIR}paw_mpelib${MODE}.o
 ${OBJDIR}paw_thermostat.o:    ${OBJDIR}paw_mpelib${MODE}.o 

 ${OBJDIR}paw_ldaplusu.o:      ${OBJDIR}paw_periodictable.o
 ${OBJDIR}paw_ioroutines.o:    ${OBJDIR}paw_periodictable.o

 ${OBJDIR}paw_ioroutines.o:    ${OBJDIR}paw_linkedlist.o
 ${OBJDIR}paw_constraints.o:   ${OBJDIR}paw_linkedlist.o
 ${OBJDIR}paw_occupations.o:   ${OBJDIR}paw_linkedlist.o
 ${OBJDIR}paw_waves.o:         ${OBJDIR}paw_linkedlist.o
 ${OBJDIR}paw_efg.o:           ${OBJDIR}paw_linkedlist.o

 ${OBJDIR}paw_atoms.o:         ${OBJDIR}paw_report.o
 ${OBJDIR}paw_ionew.o:         ${OBJDIR}paw_report.o
 ${OBJDIR}paw_occupations.o:   ${OBJDIR}paw_report.o
 ${OBJDIR}paw_qmmm.o:          ${OBJDIR}paw_report.o
 ${OBJDIR}paw_classical.o:     ${OBJDIR}paw_report.o
 ${OBJDIR}paw_thermostat.o:    ${OBJDIR}paw_report.o

 ${OBJDIR}paw_ioroutines.o:    ${OBJDIR}paw_clock.o 
 ${OBJDIR}paw_ionew.o:         ${OBJDIR}paw_waves.o

 ${OBJDIR}paw.o:               ${OBJDIR}paw_continuum.o
endif


