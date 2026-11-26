.PHONY: libxc-device-test
libxc-device-test: libpawfull.a
	mkdir -p unit-tests
	etc/f90pp $(CPPFLAGS) < $(BASEDIR)/tests/unittests/libxc_device/grid.f90 > unit-tests/libxc_grid.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/libxc_grid.x unit-tests/libxc_grid.f90 libpawfull.a $(LIBS)
	CPPAW_GPU_MODE=resident ./unit-tests/libxc_grid.x VWN
	CPPAW_GPU_MODE=resident ./unit-tests/libxc_grid.x PBE
