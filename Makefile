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

STRIP = true #strip
MODFLAGS= $(MOD_FLAG)$(TOPDIR)/Modules \
		  $(MOD_FLAG)$(TOPDIR)/LAXlib \
		  $(MOD_FLAG)$(TOPDIR)/UtilXlib \
          $(MOD_FLAG)$(TOPDIR)/FFTXlib/src \
		  $(MOD_FLAG)$(TOPDIR)/PW/src \
		  $(MOD_FLAG)$(TOPDIR)/LR_Modules \
		  $(MOD_FLAG)$(TOPDIR)/PHonon/PH \
	      $(MOD_FLAG)$(TOPDIR)/cheby/src 

IFLAGS+=-I$(TOPDIR)/FFTXlib -I$(TOPDIR)/include -I$(TOPDIR)/PHonon/PH 
DFLAGS+= -D__HASVTRIG

PHOBJS = $(TOPDIR)/PHonon/PH/libph.a $(TOPDIR)/PHonon/PH/libphaux.a
LRMODS = $(TOPDIR)/LR_Modules/liblrmod.a
PWOBJS = $(TOPDIR)/PW/src/libpw.a
QEMODS = $(TOPDIR)/Modules/libqemod.a $(TOPDIR)/KS_Solvers/Davidson/libdavid.a $(TOPDIR)/KS_Solvers/CG/libcg.a \
         $(TOPDIR)/FFTXlib/src/libqefft.a $(TOPDIR)/LAXlib/libqela.a $(TOPDIR)/UtilXlib/libutil.a $(TOPDIR)/upflib/libupf.a

EXTOBJS = $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)

# ifdef SCALAPACK_LIBS
# MINPACK = lapackified-libs distributed-libs
# EXTOBJS = ../src/libd3q.a ../minpack/distributed/libminpack.a ../minpack/lapackified/libminpack.a \
# 	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
# else
# MINPACK = lapackified-libs
# EXTOBJS = ../src/libd3q.a ../minpack/lapackified/libminpack.a \
# 	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
# endif

TLDEPS= bindir mods libs pw-lib lrmods ph-lib

tldeps:
	test -n "$(TLDEPS)" && ( cd $(TOPDIR) ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

cheby-libs: tldeps
	test -n src && ( cd src ; $(MAKE) $(MFLAGS) libcheby.a || exit 1) || :

cheby-test-libs: cheby-libs
	test -n test && ( cd test ; $(MAKE) $(MFLAGS) libchebytest.a|| exit 1) || :

checko: cheby-test-libs
	(cd test && make PROGRAM_check.o)

# d3q-libs :
# 	( cd ../D3Q/src/ && make libd3q.a)

thermal2-libs : 
	(cd $(TOPDIR)/D3Q/thermal2 && make libthermal2.a)

# lapackified-libs :
# 	( cd ../D3Q/minpack/lapackified/ && make)

# distributed-libs :
# 	( cd ../D3Q/minpack/distributed/ && make)


check : checko tldeps thermal2-libs cheby-test-libs
	$(LD) $(LDFLAGS) -o ./test/check.x test/PROGRAM_check.o test/libchebytest.a src/libcheby.a $(TOPDIR)/D3Q/thermal2/libthermal2.a  $(EXTOBJS) $(QELIBS)\
	&& ./test/check.x

# test: cheby-test-libs
# 	test -n test && (cd test; $(MAKE) $(MFLAGS) check || exit 1) || :


