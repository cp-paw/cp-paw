#!/bin/bash
################################################################################
##  paw_srclist
################################################################################
#-------------------------------------------------------------------------------
#  help message
#-------------------------------------------------------------------------------
export USAGE="\n"
USAGE="$USAGE Usage of $0:\n"
USAGE="$USAGE \t paw_srclist.sh options \n"
USAGE="$USAGE writes lists sources for the cppaw installation\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options \n"
USAGE="$USAGE \t -a list of administration files\n"
USAGE="$USAGE \t -t list of tools (w/o scripts) \n"
USAGE="$USAGE \t -p list of sources for paw.x\n"
USAGE="$USAGE \t -l list of routines for libpaw.a\n"
USAGE="$USAGE \t -v verbose (false)\n"
USAGE="$USAGE \t -h prints this help message\n"
USAGE="$USAGE \n"

#-------------------------------------------------------------------------------
#-- object files derived from FLIBLIST and  will become part                  --
#-- of the paw library libpaw.a.                                              --
#-------------------------------------------------------------------------------
export LIBLIST=" \
        paw_trace \
        paw_error \
        paw_filehandler \
        paw_clock \
        paw_timing \
        paw_linkedlist \
        paw_strings \
        paw_classical \
        paw_forcefield \
        paw_dft \
        paw_dftaddendum \
        paw_constants \
        paw_spherical \
        paw_generalpurpose \
        paw_report \
        paw_periodictable \
        paw_radial \
        paw_schroedinger \
        paw_atomlib \
        paw_specialfunctions \
        paw_usage \
        paw_selftest \
        paw_strcio \
        paw_cell \
        paw_pdos \
        paw_library \
        paw_polynom \
        paw_lmtobasics \
        paw_brillouin \
        paw_gaussian \
        paw_mpelib \
        paw_ci \
        paw_version"

#-------------------------------------------------------------------------------
#--  objects specific for the simulation code                                 --
#-------------------------------------------------------------------------------
export PAWLIST=" \
	paw_driver \
        paw_thermostat \
        paw_isolate \
        paw_dimer  \
        paw_debug \
        paw_assist \
        paw_lists \
        paw_constraints \
        paw_qmmm \
        paw_lmto \
        paw_simplelmto \
        paw_banddata \
        paw_dmft \
        paw_fft \
        paw_augmentation \
        paw_softcore \
        paw_atoms \
        paw \
        paw_efg \
        paw_ioroutines \
        paw_iotra \
        paw_ionew \
        paw_cosmo \
        paw_warmup \
        paw_optfric \
        paw_setups \
        paw_potential \
        paw_occupations \
        paw_pairpotential \
        paw_graphics \
        paw_waves1 \
        paw_waves2 \
        paw_mixer \
        paw_cg \
        paw_kpoints \
        paw_vext \
        paw_vdw \
        paw_opteels" 

# TOOLS contains the names of the source files for the tools including 
# the path relative to the tool directory $(BASEDIR)/src/
export TOOLLIST=" \
        Tools/FromPOSCAR/paw_fromposcar \
        Tools/Grab/paw_grab \
        Tools/Murnaghan/paw_murnaghan \
        Tools/PDoS/paw_dos \
        Tools/PDoS/paw_dosplot \
        Tools/Polyhedra/paw_polyhedra \
        Tools/Preopt/paw_preopt \
        Tools/Stpa/paw_stpa \
        Tools/Stpa/paw_stpreport \
        Tools/Strc/paw_strc \
        Tools/Strc/paw_tostrc \
        Tools/Strc/paw_toxyz \
        Tools/Tra/paw_cleantra \
        Tools/Tra/paw_converttra \
        Tools/Tra/paw_tra \
        Tools/Wave/paw_1davpot \
        Tools/Wave/paw_cmcwave \
        Tools/Wave/paw_wave \
        Tools/Bands/paw_bands"

# files from $(BASEDIR)/src/Buildtools to be placed in etc
export ADMINLIST=" \
        Buildtools/f90pp.in \
        Buildtools/f90pp.sed \
        Buildtools/paw_dollar_ok.sh \
        Buildtools/f90pp_tmplts.f90 \
        Buildtools/parmfilewriter.f90 \
        Buildtools/paw_versioninfo.sh" 

#-------------------------------------------------------------------------------
#--  strip whitespace                                                         --
#-------------------------------------------------------------------------------
ADMINLIST=$(echo ${ADMINLIST} | tr -s '[:blank:]') # strip extra spaces
LIBLIST=$(echo ${LIBLIST} | tr -s '[:blank:]') # strip extra spaces
PAWLIST=$(echo ${PAWLIST} | tr -s '[:blank:]') # strip extra spaces
TOOLLIST=$(echo ${TOOLLIST} | tr -s '[:blank:]') # strip extra spaces

#-------------------------------------------------------------------------------
#  Resolve arguments
#-------------------------------------------------------------------------------
while getopts :atplh OPT ; do
  case $OPT in
    l) # sources on libpaw.a
       echo "$LIBLIST" 
       ;;
    p) #sources on paw.x 
       echo "$PAWLIST"
       ;;
    t) #sources for tools
       echo "$TOOLLIST"
       ;;
    a) #sources for tools
       echo "$ADMINLIST"
       ;;
    h) echo -e $USAGE ; exit 0  ;;
    \?)   # unknown option (placed into OPTARG, if OPTSTRING starts with :)
      echo "error in $0" >&2
      echo "invalid option -$OPTARG" >&2
      echo "retrieve argument list with:" >&2
      echo "$0 -h" >&2
      exit 1
      ;;
    :)    # no argument passed
      ;;
  esac
done

