.PHONY: force-multiplier-test
force-multiplier-test: paw.x
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/force_multiplier.x $(BASEDIR)/tests/unittests/force_multiplier/force_multiplier.f90 $(ARGLIST) $(LIBS)
	./unit-tests/force_multiplier.x
