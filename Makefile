# Generated automatically from Makefile.in by configure.
#
#
SHELL=/bin/sh
#
#########################################################################
##                                                                     ##
##  make an paw executable and a paw library using gmake               ##
##  ====================================================               ##
##                                                                     ##
#########################################################################
#########################################################################
#
# The following block of variables is set by the configure
# script (AC_SUBST @...@ and sed-Commands var...)
#
#___________________________ Fortran90 compiler + compiler flags 
COMPILE=f90 -c -YEXT_NAMES=LCS -YEXT_SFX=_ 
#___________________________ Fortran77 compiler + compiler flags (for blas) 
COMPILE77=f90 -c -YEXT_NAMES=LCS -YEXT_SFX=_
#___________________________ root directory of the PAW distribution
PAWDIR=/home/ptpb/Tree/PAW
#___________________________ path to the object files
OBJDIR=/home/ptpb/Tree/PAW/bin/intel/Objects
#___________________________ path to the source files
SRCDIR=/home/ptpb/Tree/PAW/src
#___________________________ path to the binaries
BINDIR=/home/ptpb/Tree/PAW/bin
#___________________________ the PAW-binary (inclusive path)
PAWX=/home/ptpb/Tree/PAW/bin/intel
#___________________________ set to yes, if FFTW is used
FFTW=yes
#___________________________ root directory of FFTW library
FFTWDIR=/home/ptpb/Tree/PAW/lib/fftw-intel
#___________________________ root directory of MPICH library
MPICHDIR=/home/ptpb/Tree/PAW/lib/mpich-intel
#___________________________ the extension of files to be compiled
FEXT=f90
#___________________________ the hardware architecture
ARCH=intel
#___________________________ set to T, if this is Makefile.parallel 
PARALLEL=varPARALLEL
#___________________________ controls removal of preprocessed sources
OPTRM=varOPTRM
#___________________________ path to f90pp sources
F90PPDIR=${SRCDIR}/F90PP
#___________________________ directory name for not optimized objects
NONEDIR=none
#___________________________ directory name for optimized objects
FASTDIR=fast
#___________________________ directory name for parallel objects (opt.)
PARALLELDIR=fast
#___________________________ extension for parallel Makefile
PARALLELNAME=parallel
#___________________________ directory name objects of the debug-binary
DEBUGDIR=debug
#___________________________ path to the tool-directory
TOOLDIR=/home/ptpb/Tree/PAW/src/Tools
#___________________________ the way to compile the tools 
TOOLCOMP=fast
#___________________________ the GNU make utility
MAKE=gmake
#
##########################################################################
##  define parameters                                                   ##
##########################################################################
#
# at this point a distinction between parallel and sequential Makefile is
# made, as differen preprocessor variables have to be used and different
# libraries are linked
#
ifeq (${PARALLEL},T)   #select parallel versus scalar paw_mpelib.f90
  MODE = _p
  LLIBS=-lU77 -lm  -lfftw -llapack -lf77blas -latlas  -lmpich -lfmpich 
  LINK=f90   -I${OBJDIR} -L${OBJDIR}  -L/home/ptpb/Tree/PAW/lib/fftw-intel/fftw/.libs/ -L/home/ptpb/Tree/PAW/lib/ATLAS/lib/Linux_PIII -L/home/ptpb/Tree/PAW/lib/mpich-intel/lib
  F90PPOPT=-DCPPVARIABLE_PARALLEL  -DCPPVAR_U77 -DCPPVAR_FFT_FFTW -arch=${ARCH}
  #
  # MPIF90.H is an include file for the MPICH-library
  # for parallel compilation, this file has to be copied into
  # the object directory
  #
  ifneq (${ARCH},ibm)
    ifneq (${MPICHDIR},notset)
      MPICHINCLUDE=cp ${MPICHDIR}/include/mpif.h ${OBJDIR}/mpif90.h
    endif
  endif
else
  MODE = _s
  LLIBS=-lU77 -lm  -lfftw -llapack -lf77blas -latlas  
  LINK=f90   -I${OBJDIR} -L${OBJDIR}  -L/home/ptpb/Tree/PAW/lib/fftw-intel/fftw/.libs/ -L/home/ptpb/Tree/PAW/lib/ATLAS/lib/Linux_PIII 
  F90PPOPT= -DCPPVAR_U77 -DCPPVAR_FFT_FFTW -arch=${ARCH}
endif
#
# out of the file src/paw_mpelib.f a parallel _p and
# a sequential _s version are produced
#
MPE = ${OBJDIR}/paw_mpelib$(MODE).o
#
# FFTW_F77.I is an include file for the FFTW-library
# in case FFTW is linked, this file has to be copied into the 
# object directory
# 
ifeq (${FFTW},yes)
  FFTWINCLUDE=cp ${FFTWDIR}/fortran/fftw_f77.i ${OBJDIR}/FFTW_F77.I
endif
#
# paw_mpelib.f holds the subroutines necessary for parallel execution
# of the code. Via preprocessor variables, this source is split up into
# two files, one for parallel and one for sequential execution. ($MODE is
# defined in the if-statement above)
# libpaw.a holds the object files and can be linked to any code which 
# needs a PAW subroutine (leave the # after libpaw.a)
# MPE   = ${OBJDIR}/paw_mpelib$(MODE).o
PAWLIB= ${OBJDIR}/libpaw.a#         
#
#########################################################################
## Definition of the targets                                           ##
#########################################################################
#
# MIND: The Absoft and fort compiler always put the *.mod-files into the
#       current directory. In order to have them in the object-directories
#       the user types for example 'make none', wich causes a change
#       into the object-directory and second call of the Makefile.none
#       The default is then the compilation of the right binary.
# 
#
default: ${PAWX}
#
# THE PAW BINARIES
#
all: debug none fast parallel tools
#
# Two targets aim at compiling every PAW-binary (not tool). This
# is necessary, because the dependencies due to modules just have
# to be resolved at the very first calculation. In this case the
# variable ALLDEP is set to T.
#
# as the parallel objects are - by default - compiled in the
# same directory as the optimized sequential ones, there is no 
# need to resolve the dependencies. Therefore there is no
# parallel_new in the next list
#
all_new: debug_new none_new fast_new parallel tools
#
tools: 	 atom tra wave grab pdos converttra cleantra strc tostrc
#
none:   
	ALLDEP=F; export ALLDEP; cd ${OBJDIR}/${NONEDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${NONEDIR}
#
none_new:
	ALLDEP=T; export ALLDEP; cd ${OBJDIR}/${NONEDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${NONEDIR}
#
fast:   
	ALLDEP=F; export ALLDEP; cd ${OBJDIR}/${FASTDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${FASTDIR}
#
fast_new:   
	ALLDEP=T; export ALLDEP; cd ${OBJDIR}/${FASTDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${FASTDIR}
#
debug:  
	ALLDEP=F; export ALLDEP; cd ${OBJDIR}/${DEBUGDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${DEBUGDIR}
#
debug_new:  
	ALLDEP=T; export ALLDEP; cd ${OBJDIR}/${DEBUGDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${DEBUGDIR}
#
parallel:
	ALLDEP=F; export ALLDEP; cd ${OBJDIR}/${PARALLELDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${PARALLELNAME}
#
parallel_new:
	ALLDEP=T; export ALLDEP=T; cd ${OBJDIR}/${PARALLELDIR}; \
	${MAKE} -f ${PAWDIR}/Makefile.${PARALLELNAME}
#
# THE PAW TOOLS
#
atom:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_atom.x 
#
tra:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_tra.x 
#
wave:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_wave.x 
#
grab:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_grab.x 
#
pdos:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_dos.x 
#
converttra:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_converttra.x 
#
cleantra:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_cleantra.x 
#
strc:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_strc.x 
#
tostrc:
	cd ${OBJDIR}/${TOOLCOMP}; ${MAKE} -f ${PAWDIR}/Makefile.${TOOLCOMP} \
           ${BINDIR}/${ARCH}/paw_tostrc.x 
#
##########################################################################
##  define object files kept in the library                             ##
##########################################################################
#
LOBJS = \
        ${OBJDIR}/paw_trace.o \
        ${OBJDIR}/paw_error.o \
        ${OBJDIR}/paw_filehandler.o \
        ${OBJDIR}/paw_clock.o \
        ${OBJDIR}/paw_lock.o \
        ${OBJDIR}/paw_timing.o \
        ${OBJDIR}/paw_linkedlist.o \
        ${OBJDIR}/paw_strings.o \
        ${OBJDIR}/paw_dft.o \
        ${OBJDIR}/paw_dftaddendum.o \
        ${OBJDIR}/paw_constants.o \
        ${OBJDIR}/paw_spherical.o \
        ${OBJDIR}/paw_generalpurpose.o \
        ${OBJDIR}/paw_report.o \
        ${OBJDIR}/paw_periodictable.o \
        ${OBJDIR}/paw_radial.o \
        ${OBJDIR}/paw_usage.o \
        ${OBJDIR}/paw_selftest.o \
        ${OBJDIR}/paw_strcio.o \
        ${OBJDIR}/paw_cell.o \
        ${OBJDIR}/paw_library.o \
        ${OBJDIR}/paw_polynom.o 
#
##########################################################################
##  define other objects                                                ##
##########################################################################
#
OBJECTS1 = \
           ${OBJDIR}/paw_thermostat.o \
           ${OBJDIR}/paw_isolate.o \
           ${OBJDIR}/paw_assist.o \
           ${OBJDIR}/paw_ldaplusu.o \
           ${OBJDIR}/paw_lists.o \
           ${OBJDIR}/paw_constraints.o \
           ${OBJDIR}/paw_fft.o \
#
OBJECTS2 = \
           ${OBJDIR}/paw_augmentation.o \
           ${OBJDIR}/paw_atoms.o \
           ${OBJDIR}/paw.o \
           ${OBJDIR}/paw_efg.o \
           ${OBJDIR}/paw_ioroutines.o \
           ${OBJDIR}/paw_iotra.o \
           ${OBJDIR}/paw_ionew.o \
           ${OBJDIR}/paw_qmmm.o \
           ${OBJDIR}/paw_classical.o \
           ${OBJDIR}/paw_continuum.o \
           ${OBJDIR}/paw_warmup.o \
#  the following is for the QMMM of Woo and Margl
#          ${OBJDIR}/mm_paw_modules.o \
#          ${OBJDIR}/mm_paw_core_mm.o \
#          ${OBJDIR}/mm_paw_interface.o
#          ${OBJDIR}/paw_md
#
OBJECTS3 = \
           ${OBJDIR}/paw_setups.o \
           ${OBJDIR}/paw_potential.o \
           ${OBJDIR}/paw_occupations.o \
           ${OBJDIR}/paw_pairpotential.o \
           ${OBJDIR}/paw_graphics.o \
           ${OBJDIR}/paw_waves.o \
           ${OBJDIR}/paw_pdos.o \
           ${OBJDIR}/paw_kpoints.o \
           ${OBJDIR}/paw_vext.o
#          ${OBJDIR}/paw_optic.o \
#
OBJECTS= ${OBJECTS1} ${OBJECTS2} ${OBJECTS3} 
# 
#
# Definition of objects for the compilation of the PAW-tools
#
# Objects for atomic setup tool
OBJATOM = \
          ${OBJDIR}/paw_diffgl.o \
          ${OBJDIR}/paw_atom.o
# Objects for paw_tra.x
OBJTRA =  \
          ${OBJDIR}/paw_tra.o
# Objects for paw_wave.x
OBJWAVE=  \
          ${OBJDIR}/paw_wave.o
# Objects for paw_grab.x
OBJGRAB=  \
          ${OBJDIR}/paw_grab.o
# Objects for paw_dos.x
OBJPDOS=  \
          ${OBJDIR}/paw_dos.o
# Objects for paw_converttra.x
OBJCONVERTTRA=  \
          ${OBJDIR}/paw_converttra.o
# Objects for paw_cleantra.x
OBJCLEANTRA=  \
          ${OBJDIR}/paw_cleantra.o
# Objects for paw_strc.x
OBJSTRC=  \
          ${OBJDIR}/paw_strc.o
# Objects for paw_tostrc.x
OBJTOSTRC=  \
          ${OBJDIR}/paw_tostrc.o


#########################################################################
##  create preprocessor                                                ##
#########################################################################
F90PP: ${BINDIR}/${ARCH}/f90pp_tmplts.x 
${BINDIR}/${ARCH}/f90pp_tmplts.x: ${F90PPDIR}/f90pp_tmplts.f
	${BINDIR}/f90pp -notemplates ${F90PPDIR}/f90pp_tmplts.f > \
            ${OBJDIR}/f90pp_tmplts_d.${FEXT} 
	${COMPILE} -o ${OBJDIR}/f90pp_tmplts.o ${OBJDIR}/f90pp_tmplts_d.${FEXT}
	${LINK} -o ${BINDIR}/${ARCH}/f90pp_tmplts.x ${OBJDIR}/f90pp_tmplts.o
	${OPTRM} ${OBJDIR}/f90pp_tmplts_d.${FEXT}
#
#########################################################################
##  link files and produce executable PAWX                             ##
#########################################################################
#
${PAWX} : F90PP ${OBJDIR}/FFTW_F77.I ${OBJDIR}/mpif90.h ${LOBJS} \
          ${OBJECTS} ${MPE} ${OBJDIR}/liblapack.a ${OBJDIR}/libblas.a
	${LINK} -o ${PAWX} ${OBJECTS} ${LOBJS} ${MPE} ${LLIBS}
	ar -ru $(PAWLIB) $(LOBJS) ${OBJDIR}/paw_mpelib_s.o
#
# produce executables for the Tools
#
${BINDIR}/${ARCH}/paw_atom.x: ${OBJATOM} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_atom.x ${OBJATOM} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_tra.x: ${OBJTRA} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_tra.x ${OBJTRA} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_wave.x: ${OBJWAVE} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_wave.x ${OBJWAVE} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_grab.x: ${OBJGRAB} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_grab.x ${OBJGRAB} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_dos.x: ${OBJPDOS} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_dos.x ${OBJPDOS} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_converttra.x: ${OBJCONVERTTRA} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_converttra.x ${OBJCONVERTTRA} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_cleantra.x: ${OBJCLEANTRA} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_cleantra.x ${OBJCLEANTRA} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_strc.x: ${OBJSTRC} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_strc.x ${OBJSTRC} -lpaw ${LLIBS} 
#
${BINDIR}/${ARCH}/paw_tostrc.x: ${OBJTOSTRC} libpaw.a
	${LINK} -o ${BINDIR}/${ARCH}/paw_tostrc.x ${OBJTOSTRC} -lpaw ${LLIBS} 
#
# copy include files if necessary (controlled by if-statements above)
#
${OBJDIR}/FFTW_F77.I:
	${FFTWINCLUDE}
#
${OBJDIR}/mpif90.h:
	${MPICHINCLUDE}
#
#${OBJDIR}/mm_paw_f77.o 
#########################################################################
##  compile code files independent of parameter                        ##
##  and place the objects in OBJDIR                                    ##
#########################################################################
#
#__special rule for parallel interface because of preprocessor___________
${MPE} : ${OBJDIR}/%$(MODE).o: ${SRCDIR}/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*$(MODE)_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*$(MODE)_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*$(MODE)_d.${FEXT}
#
#__compile lapack routines________________________________________________
${OBJDIR}/liblapack.a : ${SRCDIR}/paw_lapack.f
	${COMPILE77} -o ${OBJDIR}/paw_lapack.o ${SRCDIR}/paw_lapack.f
	ar -ru ${OBJDIR}/liblapack.a ${OBJDIR}/paw_lapack.o
	${OPTRM} ${OBJDIR}/paw_lapack.o
#
#__compile blas routines________________________________________________
${OBJDIR}/libblas.a : ${SRCDIR}/paw_blas.f
	${COMPILE77} -o ${OBJDIR}/paw_blas.o ${SRCDIR}/paw_blas.f
	ar -ru ${OBJDIR}/libblas.a ${OBJDIR}/paw_blas.o
	${OPTRM} ${OBJDIR}/paw_blas.o
#
#__general rule for object files_________________________________________
$(OBJECTS) $(LOBJS): ${OBJDIR}/%.o: ${SRCDIR}/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
#__compile the tools____________________________________________________
${OBJATOM}: ${OBJDIR}/%.o : ${TOOLDIR}/Atom/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJTRA}: ${OBJDIR}/%.o : ${TOOLDIR}/Tra/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJWAVE}: ${OBJDIR}/%.o : ${TOOLDIR}/Wave/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJGRAB}: ${OBJDIR}/%.o : ${TOOLDIR}/Grab/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJPDOS}: ${OBJDIR}/%.o : ${TOOLDIR}/PDoS/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJCONVERTTRA}: ${OBJDIR}/%.o : ${TOOLDIR}/Tra/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJCLEANTRA}: ${OBJDIR}/%.o : ${TOOLDIR}/Tra/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJSTRC}: ${OBJDIR}/%.o : ${TOOLDIR}/Strc/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#
${OBJTOSTRC}: ${OBJDIR}/%.o : ${TOOLDIR}/Strc/%.f
	${BINDIR}/f90pp ${F90PPOPT} $< >${OBJDIR}/$*_d.${FEXT}
	${COMPILE} -o $@ ${OBJDIR}/$*_d.${FEXT}
	${OPTRM} ${OBJDIR}/$*_d.${FEXT}
#########################################################################
##  include dependencies through mod files                             ##
#########################################################################
ifeq ($(ALLDEP),T)
 ${OBJDIR}/paw_filehandler.o:   ${OBJDIR}/paw_strings.o
 ${OBJDIR}/paw_linkedlist.o:    ${OBJDIR}/paw_strings.o
 ${OBJDIR}/paw_periodictable.o: ${OBJDIR}/paw_strings.o
 ${OBJDIR}/paw_ioroutines.o:    ${OBJDIR}/paw_strings.o
#
 ${OBJDIR}/paw_timing.o:        ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_linkedlist.o:    ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_report.o:        ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_waves.o:         ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_augmentation.o:  ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw.o:               ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_constraints.o:   ${OBJDIR}/paw_mpelib${MODE}.o 
 ${OBJDIR}/paw_potential.o:     ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_pairpotential.o: ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_isolate.o:       ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_qmmm.o:          ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_thermostat.o:    ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_occupations.o:   ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_ionew.o:         ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_atoms.o:         ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_classical.o:     ${OBJDIR}/paw_mpelib${MODE}.o
 ${OBJDIR}/paw_thermostat.o:    ${OBJDIR}/paw_mpelib${MODE}.o 

 ${OBJDIR}/paw_ldaplusu.o:      ${OBJDIR}/paw_periodictable.o
 ${OBJDIR}/paw_ioroutines.o:    ${OBJDIR}/paw_periodictable.o

 ${OBJDIR}/paw_ioroutines.o:    ${OBJDIR}/paw_linkedlist.o
 ${OBJDIR}/paw_constraints.o:   ${OBJDIR}/paw_linkedlist.o
 ${OBJDIR}/paw_occupations.o:   ${OBJDIR}/paw_linkedlist.o
 ${OBJDIR}/paw_waves.o:         ${OBJDIR}/paw_linkedlist.o
 ${OBJDIR}/paw_efg.o:           ${OBJDIR}/paw_linkedlist.o

 ${OBJDIR}/paw_atoms.o:         ${OBJDIR}/paw_report.o
 ${OBJDIR}/paw_ionew.o:         ${OBJDIR}/paw_report.o
 ${OBJDIR}/paw_occupations.o:   ${OBJDIR}/paw_report.o
 ${OBJDIR}/paw_qmmm.o:          ${OBJDIR}/paw_report.o
 ${OBJDIR}/paw_classical.o:     ${OBJDIR}/paw_report.o
 ${OBJDIR}/paw_thermostat.o:    ${OBJDIR}/paw_report.o

 ${OBJDIR}/paw_ioroutines.o:    ${OBJDIR}/paw_clock.o 
 ${OBJDIR}/paw_ionew.o:         ${OBJDIR}/paw_waves.o

 ${OBJDIR}/paw.o:               ${OBJDIR}/paw_continuum.o
endif

