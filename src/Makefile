# Attempt to find make.inc in the parent directory or the grandparent directory
MAKE_INC_PATH := $(wildcard ../make.inc ../../make.inc ../../../make.inc ../../../../make.inc)

# Check if make.inc was found
ifeq ($(MAKE_INC_PATH),)
  $(error make.inc not found in .. or ../..)
endif

# Convert to absolute path and get the directory part
TOPDIR := $(dir $(abspath $(MAKE_INC_PATH)))

# Optional: include make.inc if found
include $(MAKE_INC_PATH)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)$(TOPDIR)/PW/src \
          $(MOD_FLAG)$(TOPDIR)/PHonon/PH \
          $(MOD_FLAG)$(TOPDIR)/LR_Modules \
          $(MOD_FLAG)$(TOPDIR)/D3Q/thermal2 \

SRC = $(shell ls *.f90 |grep -v PROGRAM_)
OBJ = $(shell ls *.f90 |grep -v PROGRAM_|sed -e 's/\.f90$$/.o/')

libcheby.a : $(OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@ 



