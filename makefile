# SurfSmooth3D: independent native and MATLAB builds of the imported algorithms.
# Override settings in make.inc, on the command line, or with PROFILE=<name>.
.DEFAULT_GOAL := help
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
PLATFORM := $(shell uname -s)
PROFILE ?= $(if $(filter Darwin,$(PLATFORM)),apple-silicon,linux)
ifeq ($(origin FC),default)
FC = gfortran
endif
FC ?= gfortran
CC ?= cc
AR ?= ar
NM ?= nm
FFLAGS ?= -O3 -fPIC -std=legacy
CFLAGS ?= -O3 -fPIC
OMP ?= ON
OMPFLAGS ?= -fopenmp
NATIVE_OMP_LIBS ?= -lgomp
MATLAB_OMP_LIBS ?= -lgomp
NATIVE_BLAS_LIBS ?= -lblas -llapack
MATLAB_BLAS_LIBS ?= -lmwblas -lmwlapack
EXTRA_LIBS ?= -lm -lstdc++
SHARED_EXT ?= so
SHARED_FLAGS ?= -shared -Wl,-soname,libsurfsmooth3d.so
SHARED_EXPORT_FLAGS = -Wl,--version-script,$(NATIVE)/exports.map
MATLAB_ROOT ?= /usr/local/MATLAB/R2026a
MATLAB_ARCH ?= glnxa64
MEX_EXT ?= mexa64
MEX ?= $(MATLAB_ROOT)/bin/mex
MATLAB_RUNTIME_PIN ?= 0
MEX_FLAGS ?= -compatibleArrayDims -DMWF77_UNDERSCORE1
MWRAP ?= mwrap
MWFLAGS ?= -c99complex -i8
FMM_PREFIX ?= $(ROOT)/build/deps/fmm3d-reference/install
FMM_STATIC ?= $(firstword $(wildcard $(FMM_PREFIX)/lib/libfmm3d.a $(FMM_PREFIX)/lib-static/libfmm3d.a $(FMM_PREFIX)/libfmm3d.a))
PREFIX ?= $(ROOT)/install
BUILD ?= $(ROOT)/build
STEP_BUILD ?= $(BUILD)/step
STEP_DEPENDENCY_PREFIX ?= $(ROOT)/build/deps/cad/install
STEP_MEX ?= ON
STEP_CMAKE_FLAGS ?=
-include $(ROOT)/config/$(PROFILE).mk
-include $(ROOT)/make.inc

NATIVE := $(BUILD)/native
MATLAB := $(BUILD)/matlab
ifeq ($(OMP),OFF)
FOMP :=
NATIVE_OMP_LIBS :=
MATLAB_OMP_LIBS :=
else
FOMP := $(OMPFLAGS)
endif
CORE_NAMES := ModType_Smooth_Surface Mod_TreeLRD Mod_Fast_Sigma Mod_Plot_Tools_sigma Mod_Feval Mod_Smooth_Surface Mod_Adaptive_Smooth_Surface Mod_Adaptive_Blend_Surface cisurf_loadmsh cisurf_skeleton cisurf_plottools cisurf_tritools surface_smoother tfmm_setsub
SUPPORT_NAMES := koornexps setops sort
NATIVE_OBJS := $(addprefix $(NATIVE)/obj/,$(addsuffix .o,$(CORE_NAMES) $(SUPPORT_NAMES) lapack_wrap))
MATLAB_OBJS := $(addprefix $(MATLAB)/obj/,$(addsuffix .o,$(CORE_NAMES) $(SUPPORT_NAMES) lapack_wrap_64))
NATIVE_CORE_OBJS := $(addprefix $(NATIVE)/obj/,$(addsuffix .o,$(CORE_NAMES)))
NATIVE_LIB := $(NATIVE)/lib/libsurfsmooth3d.a
NATIVE_SHARED := $(NATIVE)/lib/libsurfsmooth3d.$(SHARED_EXT)
MATLAB_LIB := $(MATLAB)/lib/libsurfsmooth3d_matlab.a
NATIVE_LINK := $(NATIVE_LIB) $(FMM_STATIC) $(NATIVE_BLAS_LIBS) $(NATIVE_OMP_LIBS) $(EXTRA_LIBS)
MEX_NAMES := surfsmooth3d_routs surfsmooth3d_adaptive_routs surfsmooth3d_adaptive_blend_routs
MEX_FILES := $(addprefix $(ROOT)/matlab/,$(addsuffix .$(MEX_EXT),$(MEX_NAMES)))
TEST_NAMES := test_sigma_modes test_newton_radius_guard test_adaptive_smoother test_adaptive_blend test_fixed_sources test_surfsmooth
TEST_BINS := $(addprefix $(NATIVE)/bin/,$(TEST_NAMES))
FIXTURES := $(ROOT)/tests/fixtures/cylinder
TEST_OUTPUT := $(NATIVE)/test-output
TEST_LOGS := $(NATIVE)/test-logs

.PHONY: help lib matlab generate-matlab test-native test-c test-matlab install install-matlab install-docs step install-step check-fmm clean FORCE
help:
	@echo "make lib | matlab | step | test-native | test-c | test-matlab"
	@echo "make install | install-matlab | install-step PREFIX=/path/to/install"
	@echo "Set FMM_PREFIX to an external FMM3D installation (lib/libfmm3d.a)."
	@echo "PROFILE=$(PROFILE), FC=$(FC), FMM_PREFIX=$(FMM_PREFIX)"

check-fmm:
	@test -n "$(FMM_STATIC)" && test -f "$(FMM_STATIC)" || { echo 'FMM3D archive not found: set FMM_PREFIX or FMM_STATIC.' >&2; exit 1; }

lib: check-fmm $(NATIVE_LIB) $(NATIVE_SHARED) $(NATIVE)/bin/smooth_surface $(NATIVE)/bin/smooth_surface_shared

# Record effective settings without changing a stamp when its contents match.
# In particular, selecting a different (possibly older) external archive or
# changing the local MATLAB runtime option must rebuild the affected outputs.
shquote = '$(subst ','"'"',$(1))'
FORCE:
$(NATIVE)/compile.settings $(MATLAB)/compile.settings: FORCE
	@mkdir -p $(@D)
	@printf '%s\n' $(call shquote,FC=$(FC)) $(call shquote,FFLAGS=$(FFLAGS)) $(call shquote,FOMP=$(FOMP)) > $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@
	@rm -f $@.tmp

$(NATIVE)/link.settings: FORCE | check-fmm
	@mkdir -p $(@D)
	@printf '%s\n' $(call shquote,FC=$(FC)) $(call shquote,FMM_STATIC=$(FMM_STATIC)) $(call shquote,NATIVE_BLAS_LIBS=$(NATIVE_BLAS_LIBS)) $(call shquote,NATIVE_OMP_LIBS=$(NATIVE_OMP_LIBS)) $(call shquote,EXTRA_LIBS=$(EXTRA_LIBS)) $(call shquote,SHARED_FLAGS=$(SHARED_FLAGS)) $(call shquote,SHARED_EXPORT_FLAGS=$(SHARED_EXPORT_FLAGS)) > $@.tmp
	@cksum $(call shquote,$(FMM_STATIC)) >> $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@
	@rm -f $@.tmp

$(MATLAB)/link.settings: FORCE | check-fmm
	@mkdir -p $(@D)
	@printf '%s\n' $(call shquote,MEX=$(MEX)) $(call shquote,FMM_STATIC=$(FMM_STATIC)) $(call shquote,MEX_FLAGS=$(MEX_FLAGS)) $(call shquote,MEX_PLATFORM_FLAGS=$(MEX_PLATFORM_FLAGS)) $(call shquote,MATLAB_BLAS_LIBS=$(MATLAB_BLAS_LIBS)) $(call shquote,MATLAB_OMP_LIBS=$(MATLAB_OMP_LIBS)) $(call shquote,MATLAB_RUNTIME_PIN=$(MATLAB_RUNTIME_PIN)) > $@.tmp
	@cksum $(call shquote,$(FMM_STATIC)) >> $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@
	@rm -f $@.tmp

$(NATIVE_OBJS): $(NATIVE)/compile.settings
$(MATLAB_OBJS): $(MATLAB)/compile.settings

# FMM3D provides legeexps and prini, exactly as in the tested original build.
# Their provenance copies are not compiled; there is one provider per link.
define profile_rules
$(1)/obj/%.o: $(ROOT)/src/multiscale_mesher/%.f90
	@mkdir -p $(1)/obj $(1)/mod
	$(FC) $(FFLAGS) $(FOMP) -J$(1)/mod -I$(1)/mod -c $$< -o $$@
$(1)/obj/%.o: $(ROOT)/src/multiscale_mesher/%.f
	@mkdir -p $(1)/obj $(1)/mod
	$(FC) $(FFLAGS) $(FOMP) -J$(1)/mod -I$(1)/mod -c $$< -o $$@
$(1)/obj/%.o: $(ROOT)/src/support/%.f90
	@mkdir -p $(1)/obj $(1)/mod
	$(FC) $(FFLAGS) $(FOMP) -J$(1)/mod -I$(1)/mod -c $$< -o $$@
$(1)/obj/%.o: $(ROOT)/src/support/%.f
	@mkdir -p $(1)/obj $(1)/mod
	$(FC) $(FFLAGS) $(FOMP) -J$(1)/mod -I$(1)/mod -c $$< -o $$@
$(1)/obj/koornexps.o: $(ROOT)/src/special_functions/koornexps.f90 $(ROOT)/src/special_functions/koorn-uvs-dat.txt $(ROOT)/src/special_functions/koorn-wts-dat.txt
	@mkdir -p $(1)/obj $(1)/mod
	$(FC) $(FFLAGS) $(FOMP) -J$(1)/mod -I$(1)/mod -I$(ROOT)/src/special_functions -c $$< -o $$@
$(1)/obj/Mod_Fast_Sigma.o: $(1)/obj/Mod_TreeLRD.o $(1)/obj/ModType_Smooth_Surface.o
$(1)/obj/Mod_Plot_Tools_sigma.o: $(1)/obj/Mod_Fast_Sigma.o
$(1)/obj/Mod_Feval.o: $(1)/obj/Mod_Fast_Sigma.o
$(1)/obj/Mod_Smooth_Surface.o: $(1)/obj/Mod_Feval.o
$(1)/obj/Mod_Adaptive_Smooth_Surface.o: $(1)/obj/Mod_Smooth_Surface.o
$(1)/obj/Mod_Adaptive_Blend_Surface.o: $(1)/obj/Mod_Adaptive_Smooth_Surface.o
$(1)/obj/surface_smoother.o: $(1)/obj/Mod_Smooth_Surface.o $(1)/obj/Mod_Plot_Tools_sigma.o
$(1)/obj/cisurf_loadmsh.o $(1)/obj/cisurf_skeleton.o: $(1)/obj/ModType_Smooth_Surface.o
$(1)/obj/cisurf_plottools.o: $(1)/obj/Mod_Smooth_Surface.o
endef
$(eval $(call profile_rules,$(NATIVE)))
$(eval $(call profile_rules,$(MATLAB)))

$(NATIVE_LIB): $(NATIVE_OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $(NATIVE_OBJS)

$(MATLAB_LIB): $(MATLAB_OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $(MATLAB_OBJS)

# Only smoother symbols are public. Vendored support and linked FMM internals
# are private in the shared library; the static archive is not co-link safe
# with another unprefixed copy of its support routines.
$(NATIVE)/exports.txt: $(NATIVE_CORE_OBJS)
	$(NM) -g $(NATIVE_CORE_OBJS) | awk '$$2 ~ /^[TDBS]$$/ { print $$3 }' | sort -u > $@
$(NATIVE)/exports.map: $(NATIVE)/exports.txt
	awk 'BEGIN {print "{ global:"} {print $$0 ";"} END {print "local: *; };"}' $< > $@
$(NATIVE_SHARED): $(NATIVE_LIB) $(NATIVE)/exports.txt $(NATIVE)/exports.map $(FMM_STATIC) $(NATIVE)/link.settings | check-fmm
	$(FC) $(SHARED_FLAGS) $(SHARED_EXPORT_FLAGS) -o $@ $(NATIVE_OBJS) $(FMM_STATIC) $(NATIVE_BLAS_LIBS) $(NATIVE_OMP_LIBS) $(EXTRA_LIBS)

$(NATIVE)/bin/smooth_surface: $(ROOT)/examples/fortran/smooth_surface.f90 $(NATIVE_LIB) $(FMM_STATIC) $(NATIVE)/link.settings | check-fmm
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) $(FOMP) -I$(NATIVE)/mod $< $(NATIVE_LINK) -o $@
$(NATIVE)/bin/smooth_surface_shared: $(ROOT)/examples/fortran/smooth_surface.f90 $(NATIVE_SHARED)
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) $(FOMP) -I$(NATIVE)/mod $< -L$(NATIVE)/lib -lsurfsmooth3d -Wl,-rpath,$(NATIVE)/lib -o $@

$(NATIVE)/bin/test_%: $(ROOT)/tests/fortran/test_%.f90 $(NATIVE_LIB) $(FMM_STATIC) $(NATIVE)/link.settings | check-fmm
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) $(FOMP) -fcheck=all -I$(NATIVE)/mod $< $(NATIVE_LINK) -o $@

$(NATIVE)/bin/test_surfsmooth_c: $(ROOT)/tests/fortran/test_surfsmooth.c $(NATIVE_SHARED) $(NATIVE)/link.settings
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -L$(NATIVE)/lib -lsurfsmooth3d -Wl,-rpath,$(NATIVE)/lib -o $@

test-c: lib $(NATIVE)/bin/test_surfsmooth_c
	@mkdir -p $(TEST_OUTPUT) $(TEST_LOGS)
	@$(NATIVE)/bin/test_surfsmooth_c $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/c-interface > $(TEST_LOGS)/c-interface.log 2>&1 || { cat $(TEST_LOGS)/c-interface.log; exit 1; }
	@$(NATIVE)/bin/smooth_surface $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/c-reference 8 0 0 > $(TEST_LOGS)/c-reference.log 2>&1 || { cat $(TEST_LOGS)/c-reference.log; exit 1; }
	@cmp $(TEST_OUTPUT)/c-interface_o08_r00.go3 $(TEST_OUTPUT)/c-reference_o08_r00.go3
	@echo "PASS: C-string and native Fortran interface outputs match."

test-native: lib $(TEST_BINS)
	@mkdir -p $(TEST_OUTPUT) $(TEST_LOGS)
	@$(NATIVE)/bin/test_sigma_modes > $(TEST_LOGS)/sigma.log 2>&1 || { cat $(TEST_LOGS)/sigma.log; exit 1; }
	@$(NATIVE)/bin/test_newton_radius_guard $(TEST_OUTPUT)/guard > $(TEST_LOGS)/newton-guard.log 2>&1 || { cat $(TEST_LOGS)/newton-guard.log; exit 1; }
	@$(NATIVE)/bin/test_adaptive_smoother $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/adaptive > $(TEST_LOGS)/adaptive.log 2>&1 || { cat $(TEST_LOGS)/adaptive.log; exit 1; }
	@$(NATIVE)/bin/test_adaptive_blend $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/blend > $(TEST_LOGS)/blend.log 2>&1 || { cat $(TEST_LOGS)/blend.log; exit 1; }
	@$(NATIVE)/bin/test_fixed_sources $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/fixed > $(TEST_LOGS)/fixed-sources.log 2>&1 || { cat $(TEST_LOGS)/fixed-sources.log; exit 1; }
	@$(NATIVE)/bin/test_surfsmooth $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/legacy > $(TEST_LOGS)/legacy-kernel.log 2>&1 || { cat $(TEST_LOGS)/legacy-kernel.log; exit 1; }
	@$(NATIVE)/bin/smooth_surface $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/static 8 0 0 > $(TEST_LOGS)/static-link.log 2>&1 || { cat $(TEST_LOGS)/static-link.log; exit 1; }
	@$(NATIVE)/bin/smooth_surface_shared $(FIXTURES)/cylinder.gidmsh $(FIXTURES)/cylinder_skeleton.txt $(TEST_OUTPUT)/shared 8 0 0 > $(TEST_LOGS)/shared-link.log 2>&1 || { cat $(TEST_LOGS)/shared-link.log; exit 1; }
	@cmp $(TEST_OUTPUT)/static_o08_r00.go3 $(TEST_OUTPUT)/shared_o08_r00.go3
	@echo "PASS: six native suites; static/shared public driver outputs match. Logs: $(TEST_LOGS)"

matlab: check-fmm $(MEX_FILES)
ifeq ($(MATLAB_RUNTIME_PIN),1)
MEX_GATEWAY = $(ROOT)/config/runtime/$(1)_gateway.c
else
MEX_GATEWAY = $(ROOT)/matlab/$(1).c
endif
define mex_rule
$(ROOT)/matlab/$(1).$(MEX_EXT): $(ROOT)/matlab/$(1).c $(MATLAB_LIB) $(FMM_STATIC) $(MATLAB)/link.settings $(wildcard $(ROOT)/config/runtime/$(1)_gateway.c) | check-fmm
	@mkdir -p $(MATLAB)/prefs
	env MATLAB_PREFDIR=$(MATLAB)/prefs $(MEX) $(call MEX_GATEWAY,$(1)) $(MATLAB_LIB) $(FMM_STATIC) $$(MEX_FLAGS) $$(MEX_PLATFORM_FLAGS) $$(MATLAB_BLAS_LIBS) $$(MATLAB_OMP_LIBS) -lm -lstdc++ -lgfortran -output $(ROOT)/matlab/$(1)
endef
$(foreach target,$(MEX_NAMES),$(eval $(call mex_rule,$(target))))

generate-matlab:
	cd $(ROOT)/matlab && $(MWRAP) $(MWFLAGS) -list -mex surfsmooth3d_routs -mb surfsmooth3d_routs.mw
	cd $(ROOT)/matlab && $(MWRAP) $(MWFLAGS) -mex surfsmooth3d_routs -c surfsmooth3d_routs.c surfsmooth3d_routs.mw

install: lib install-docs
	install -d $(PREFIX)/lib $(PREFIX)/include/surfsmooth3d $(PREFIX)/bin
	install -m 644 $(NATIVE_LIB) $(NATIVE_SHARED) $(PREFIX)/lib/
	install -m 644 $(NATIVE)/mod/*.mod $(PREFIX)/include/surfsmooth3d/
	install -m 755 $(NATIVE)/bin/smooth_surface $(PREFIX)/bin/surfsmooth3d

install-matlab: matlab install-docs
	install -d $(PREFIX)/matlab
	cp -R $(ROOT)/matlab/+surfsmooth3d $(PREFIX)/matlab/
	cp $(ROOT)/matlab/setup_surfsmooth3d.m $(MEX_FILES) $(PREFIX)/matlab/

install-docs:
	install -d $(PREFIX)/share/SurfSmooth3D
	cp $(ROOT)/README.md $(ROOT)/LICENSE $(ROOT)/COPYRIGHT $(ROOT)/THIRD_PARTY.md $(PREFIX)/share/SurfSmooth3D/
	cp -R $(ROOT)/docs $(ROOT)/provenance $(PREFIX)/share/SurfSmooth3D/

TEST_GROUP ?= all
TEST_THREADS ?= 2
test-matlab: matlab
	@mkdir -p $(BUILD)/matlab-test-prefs
	env LC_ALL=C LANG=C OMP_NUM_THREADS=$(TEST_THREADS) MATLAB_PREFDIR=$(BUILD)/matlab-test-prefs $(MATLAB_ROOT)/bin/matlab -batch "addpath('$(ROOT)/tests/matlab'); run_surfsmooth3d_tests('$(TEST_GROUP)');"

step:
	cmake -S $(ROOT)/step_mesher -B $(STEP_BUILD) -DCMAKE_INSTALL_PREFIX=$(PREFIX) -DSTEP_MESHER_DEPENDENCY_PREFIX=$(STEP_DEPENDENCY_PREFIX) -DSTEP_MESHER_BUILD_MEX=$(STEP_MEX) -DMatlab_ROOT_DIR=$(MATLAB_ROOT) $(STEP_CMAKE_FLAGS)
	cmake --build $(STEP_BUILD) --parallel

install-step:
	cmake --install $(STEP_BUILD) --prefix $(PREFIX)

clean:
	rm -rf $(NATIVE) $(MATLAB)
