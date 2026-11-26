XC_FUNCTIONAL ?= SCAN
XC_MODEL ?=
XC_MODEL_PRECISION ?= published
XC_DEVICE ?= CPU

.PHONY: xc-spin-closed-test xc-spin-open-test xc-spin-noncollinear-test

unit-tests/xc_spin.x: libpawfull.a $(BASEDIR)/tests/unittests/xc_spin/collinear.f90
	mkdir -p unit-tests
	etc/f90pp $(CPPFLAGS) < $(BASEDIR)/tests/unittests/xc_spin/collinear.f90 > unit-tests/xc_spin.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o $@ unit-tests/xc_spin.f90 libpawfull.a $(LIBS)

xc-spin-closed-test: unit-tests/xc_spin.x
	./unit-tests/xc_spin.x "$(XC_FUNCTIONAL)" "$(XC_MODEL)" "$(XC_MODEL_PRECISION)" closed

xc-spin-open-test: unit-tests/xc_spin.x
	./unit-tests/xc_spin.x "$(XC_FUNCTIONAL)" "$(XC_MODEL)" "$(XC_MODEL_PRECISION)" open

xc-spin-noncollinear-test: libpawfull.a
	mkdir -p unit-tests
	etc/f90pp $(CPPFLAGS) < $(BASEDIR)/tests/unittests/xc_spin/noncollinear.f90 > unit-tests/xc_noncollinear.f90
	$(LD) $(FCFLAGS) $(LDFLAGS) -I. -o unit-tests/xc_noncollinear.x unit-tests/xc_noncollinear.f90 libpawfull.a $(LIBS)
	./unit-tests/xc_noncollinear.x "$(XC_FUNCTIONAL)" "$(XC_MODEL)" "$(XC_MODEL_PRECISION)" "$(XC_DEVICE)"
