# Tuning parameters.
MASS_DIFF_P="0.2"
COM_CORR="com" # turn on with 'com'

# All test files
testinput := $(shell find 'testinput' -name '*.png')
testnames = $(testinput:testinput/%.png=%)

# The default file and all test files.
allexperiments = default $(testnames)

all: run/default.plots

all_base: $(allexperiments:%=run/%.run0)

all_expt: $(allexperiments:%=run/%.plots)

# Don't delete my intermediate files!
.SECONDARY:

adaptive: adaptive.f95
	mkdir -p bin
	gfortran -O5 adaptive.f95 -o adaptive

# Always run, even if there is a file called 'clean'
.PHONY: clean
clean:
	rm -f bin/*
	rm -f testinput/*.data
	rm -f run/*.run*

clean_expt:
	rm -f run/*.run run/*.plots run/*/*.data

run/default/density_project.data:
	mkdir -p run/default
	cp -f density_project.data run/default/

run/%/density_project.data: testinput/%.png
	mkdir -p run/$*
	python testinput/totxt.py testinput/$*.png > $@

run/%.run: adaptive run/%/density_project.data
	./adaptive 6 0 $(MASS_DIFF_P) $(COM_CORR) run/$*/density_project.data run/$*/ | tee $@

run/%.run0: adaptive run/%/density_project.data
	mkdir -p run/$*/level0
	./adaptive 0 0 0.0 "" run/$*/density_project.data run/$*/level0/ | tee $@

run/%.plots: run/%.run
	python graphing.py --prefix "run/$*" --save --all | tee $@
	python graphing.py --prefix "run/$*" --save --all --polar | tee $@

# Only used for development.
run/%.test: run/%.run
	python graphing.py --prefix "run/$*" -d2
