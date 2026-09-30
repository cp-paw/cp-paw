.PHONY: hartree-strain-test
hartree-strain-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/hartree_strain/hartree_strain.f90 > unit-tests/hartree_strain.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/hartree_strain.x unit-tests/hartree_strain.f90 $(ARGLIST) $(LIBS)
	./unit-tests/hartree_strain.x
