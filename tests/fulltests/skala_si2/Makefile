SHELL=/bin/sh

.PHONY: all check clean

all: check model.fun
	@echo "running Skala PAW integration test..."
	@export PAWDIR="not defined"; ${PAWX} skala_si2.cntl 2>err.log >out.log
	@${SHELL} analyse.sh

check:
	@test -n "${SKALA_MODEL}" || { \
	  echo "SKALA_MODEL must name a Skala 1.1 TorchScript model" >&2; exit 2; }
	@test -f "${SKALA_MODEL}" || { \
	  echo "SKALA_MODEL does not exist: ${SKALA_MODEL}" >&2; exit 2; }

model.fun: check
	@ln -sf "$$(realpath "${SKALA_MODEL}")" $@

clean:
	@rm -f model.fun err.log out.log skala_si2.banddata skala_si2.pdos \
	  skala_si2.prot skala_si2.rstrt skala_si2.strc_out \
	  skala_si2_constr.report skala_si2_stpforz14.myxml
