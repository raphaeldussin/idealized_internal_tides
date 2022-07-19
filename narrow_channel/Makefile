NANGLES := 24
ANGLE := $(shell seq 1 6 ${NANGLES})
ROTATION := 0 1

#JOBS := $(addprefix regr,${ANGLE})
#.PHONY: all ${JOBS}
#all: ${JOBS} ; echo "mkdir -p regression/$@ success"
#${JOBS}: regr%: ; echo "mkdir -p regression/$*"

#.PHONY: regr

#regr:
#	$(foreach a, $(ANGLE), $(foreach r, $(ROTATION), $(shell mkdir -p regressions/regr_angle$(a)_rot$(r))))


#all: regression/regr_angle$(ANGLE)/Makefile

#.PRECIOUS: $(foreach a,$(ANGLE),regression/regr_angle$(a)/Makefile)

tmp = $(addprefix regression/regr_angle,$(ANGLE))
tmp2 = $(addsuffix $(tmp), _rot)
tmp3 = $(addsuffix $(tmp2), $(ROTATION))

#all: $(addprefix regression/regr_angle,$(ANGLE))
all: ${tmp3}


#all: regression/regr_angle${ANGLE}/Makefile

regression/regr_angle%:
	@\mkdir -p regression/regr_angle$*

regression/regr_angle%_rot%:
	@\mkdir -p regression/regr_angle$*_rot$*

#regression/regr_angle%/Makefile:
#	@\mkdir -p $(@D)
#	@touch $(@D)/Makefile

#regr%:
#	echo "mkdir -p regression/$*"

clean:
	@rm -f available_diags.000000 time_stamp.out ocean.stats{,.nc} logfile.000000.out exitcode CPU_stats core

cleannc:
	@rm -f *.nc
	@rm -rf INPUT RESTART
