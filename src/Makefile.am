AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../kalign -I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  ../kalign/libkalign.la ../biomcmc-lib/lib/libbiomcmc.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = tatajuba 
tatajuba_SOURCES = tatajuba.h main.c kseq.h
tatajuba_LDADD = $(LDADD)
