.PHONY: skala-reconstruction-test
skala-reconstruction-test: libpaw.a
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o skala_reconstruction_test.x $(BASEDIR)/tests/unittests/skala_reconstruction/skala_reconstruction.f90 libpaw.a $(LIBS)
	./skala_reconstruction_test.x
