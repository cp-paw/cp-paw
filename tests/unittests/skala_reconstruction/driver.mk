.PHONY: skala-reconstruction-test
skala-reconstruction-test: libpaw.a
	mkdir -p unit-tests
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/skala_reconstruction.x $(BASEDIR)/tests/unittests/skala_reconstruction/skala_reconstruction.f90 libpaw.a $(LIBS)
	./unit-tests/skala_reconstruction.x
