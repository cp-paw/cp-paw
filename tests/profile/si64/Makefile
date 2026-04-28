SHELL=/bin/sh
#
TEST ?= si64
PROFILE_PREFIX ?= $(TEST)_accel_profile
STP ?= ../../fulltests/si2/stp.cntl
NVBLAS ?= no
ifeq ($(NVBLAS),yes)
  NVBLAS_WRAPPER ?= paw_nvblas.sh --
else
  NVBLAS_WRAPPER ?=
endif
RUN_OUTPUTS = $(PROFILE_PREFIX).csv $(PROFILE_PREFIX).rank*.csv \
              $(TEST).prot $(TEST).rstrt $(TEST).strc_out \
              $(TEST)_r.tra $(TEST).banddata $(TEST).pdos \
              $(TEST)_constr.report $(TEST)_stpforz*.myxml \
              out out.log err.log nvblas.log
#
.PHONY: all summary benchmark nsys overnight clean
#
all:
	@test -n "$${PAWX}" || { echo "Set PAWX to the CP-PAW command or MPI launcher."; exit 2; }
	@cp $(STP) stp.cntl
	@rm -f $(RUN_OUTPUTS)
	@echo "running $(TEST) profile..."
	@CPPAW_ACCEL_PROFILE_FILE=$(PROFILE_PREFIX) $(NVBLAS_WRAPPER) $${PAWX} $(TEST).cntl 2> err.log > out.log
	@$(MAKE) --no-print-directory summary

summary:
	@python3 profile_summary.py "$(PROFILE_PREFIX)"*.csv

benchmark:
	@./run_benchmark.sh

nsys:
	@./run_nsys.sh

overnight:
	@./run_overnight.sh

clean:
	@rm -f $(RUN_OUTPUTS) stp.cntl
