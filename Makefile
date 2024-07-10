# Makefile for D3Q
# Adapted from TDDFPT main Makefile

include ../../make.inc

STRIP = true #strip
MODFLAGS= $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../UtilXlib \
          $(MOD_FLAG)../../FFTXlib/src $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules $(MOD_FLAG)../../PHonon/PH \
	  $(MOD_FLAG)../src $(MOD_FLAG). $(MOD_FLAG)../thermal2

IFLAGS+=-I../../FFTXlib -I../../include -I../../PHonon/PH 
DFLAGS+= -D__HASVTRIG

PHOBJS = ../../PHonon/PH/libph.a ../../PHonon/PH/libphaux.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/Davidson/libdavid.a ../../KS_Solvers/CG/libcg.a \
         ../../FFTXlib/src/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../upflib/libupf.a

EXTOBJS = ../src/libd3q.a $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)

ifdef SCALAPACK_LIBS
MINPACK = lapackified-libs distributed-libs
EXTOBJS = ../src/libd3q.a ../minpack/distributed/libminpack.a ../minpack/lapackified/libminpack.a \
	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
else
MINPACK = lapackified-libs
EXTOBJS = ../src/libd3q.a ../minpack/lapackified/libminpack.a \
	  $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
endif

TLDEPS= bindir mods libs pw-lib lrmods ph-lib

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

cheby-libs: tldeps
	test -n src && ( cd src ; $(MAKE) $(MFLAGS) libcheby.a || exit 1) || :

cheby-test-libs: cheby-libs
	test -n test && ( cd test ; $(MAKE) $(MFLAGS) libchebytest.a|| exit 1) || :

checko: cheby-test-libs
	(cd test && make PROGRAM_check.o)

d3q-libs :
	( cd ../src/ && make libd3q.a)

thermal2-libs : 
	(cd ../thermal2 && make libthermal2.a)

lapackified-libs :
	( cd ../minpack/lapackified/ && make)

distributed-libs :
	( cd ../minpack/distributed/ && make)


check : checko tldeps d3q-libs thermal2-libs cheby-test-libs $(MINPACK)
	$(LD) $(LDFLAGS) -o ./test/check.x test/PROGRAM_check.o test/libchebytest.a src/libcheby.a ../thermal2/libthermal2.a  $(EXTOBJS) $(QELIBS)\
	&& ./test/check.x

test: cheby-test-libs
	test -n test && (cd test; $(MAKE) $(MFLAGS) check || exit 1) || :


