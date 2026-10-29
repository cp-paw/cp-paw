.PHONY: scan-test

scan-test: libpaw.a
	mkdir -p unit-tests
	etc/f90pp $(CPPFLAGS) < $(BASEDIR)/tests/unittests/scan/test_scan.f90 > unit-tests/test_scan.f90
	etc/f90pp $(CPPFLAGS) < $(BASEDIR)/tests/unittests/scan/reference_scan.f90 > unit-tests/reference_scan.f90
	etc/f90pp $(CPPFLAGS) < $(BASEDIR)/tests/unittests/scan/block_scan.f90 > unit-tests/block_scan.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/test_scan.x unit-tests/test_scan.f90 unit-tests/reference_scan.f90 unit-tests/block_scan.f90 libpaw.a $(LIBS)
	./unit-tests/test_scan.x
