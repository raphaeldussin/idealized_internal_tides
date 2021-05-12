# Include any local configurations
-include config.mk

TEMPLATE ?= src/mkmf/templates/cheyenne-intel.mk

FMS_MODS = affinity astronomy axis_utils block_control constants data_override diag_manager exchange field_manager fms horiz_interp include interpolator memutils mosaic mpp oda_tools platform random_numbers time_interp time_manager
FMS_MODS = affinity amip_interp astronomy axis_utils block_control column_diagnostics constants coupler data_override diag_integral diag_manager drifters exchange fft field_manager fms fv3gfs horiz_interp include interpolator memutils monin_obukhov mosaic mpp oda_tools platform random_numbers sat_vapor_pres station_data test_fms time_interp time_manager topography tracer_manager tridiagonal


MOM6_PATHS = config_src/drivers/solo_driver config_src/memory/dynamic_symmetric config_src/external config_src/infra/FMS1 src

build/MOM6: build/Makefile
	cd $(@D) ; make
build/Makefile: build/path_names
	cd $(@D) ; ../src/mkmf/bin/mkmf -t ../$(TEMPLATE) -p MOM6 -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
build/path_names: $(shell find src/MOM6 -name "*.F90")
	@mkdir -p $(@D)
	cd $(@D) ; ../src/mkmf/bin/list_paths -l $(foreach m,$(FMS_MODS),../src/FMS/$(m)) $(foreach m,$(MOM6_PATHS),../src/MOM6/$(m))
	cd $(@D) ; sed -i '/\/test_/d' path_names

clean:
	rm -rf build
