.PHONY: orthogonality-tolerance-test
orthogonality-tolerance-test: paw.x
	mkdir -p unit-tests
	$(F90PP) $(CPPFLAGS) < $(BASEDIR)/tests/unittests/orthogonality_tolerance/orthogonality_tolerance.f90 > unit-tests/orthogonality_tolerance.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/orthogonality_tolerance.x unit-tests/orthogonality_tolerance.f90 $(ARGLIST) $(LIBS)
	./unit-tests/orthogonality_tolerance.x
	@for mode in nan inf loose tight zero; do \
	  if GFORTRAN_UNBUFFERED_ALL=1 ./unit-tests/orthogonality_tolerance.x $$mode > unit-tests/orthogonality_$$mode.log 2>&1; then exit 1; fi; \
	  grep -q 'ORTHOTOL MUST' unit-tests/orthogonality_$$mode.log || exit 1; \
	done
