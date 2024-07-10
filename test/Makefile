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
          $(MOD_FLAG)$(TOPDIR)/cheby/src \

# PHOBJS = $(TOPDIR)PHonon/PH/libph.a $(TOPDIR)PHonon/PH/libphaux.a
# LRMODS = $(TOPDIR)LR_Modules/liblrmod.a
# PWOBJS = $(TOPDIR)PW/src/libpw.a
# QEMODS = $(TOPDIR)Modules/libqemod.a $(TOPDIR)KS_Solvers/Davidson/libdavid.a $(TOPDIR)KS_Solvers/CG/libcg.a \
#          $(TOPDIR)FFTXlib/src/libqefft.a $(TOPDIR)LAXlib/libqela.a $(TOPDIR)UtilXlib/libutil.a $(TOPDIR)upflib/libupf.a


# ifdef SCALAPACK_LIBS
# MINPACK = lapackified-libs distributed-libs
# EXTOBJS = $(TOPDIR)/D3Q/src/libd3q.a $(TOPDIR)/D3Q/minpack/distributed/libminpack.a $(TOPDIR)/D3Q/minpack/lapackified/libminpack.a \
# 	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
# else
# MINPACK = lapackified-libs
# EXTOBJS = $(TOPDIR)/D3Q/src/libd3q.a $(TOPDIR)/D3Q/minpack/lapackified/libminpack.a \
# 	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
# endif

# LIBOBJS = $(TOPDIR)/FFTXlib/src/libqefft.a $(TOPDIR)/LAXlib/libqela.a

SRC = $(shell ls *.f90 |grep -v PROGRAM_)
OBJ = $(shell ls *.f90 |grep -v PROGRAM_|sed -e 's/\.f90$$/.o/')

libchebytest.a : $(OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@ 