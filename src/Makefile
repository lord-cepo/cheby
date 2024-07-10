include ../../../make.inc


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



