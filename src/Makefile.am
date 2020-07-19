AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../kalign -I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  ../kalign/libkalign.la ../biomcmc-lib/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = tatajuba 
tatajuba_SOURCES = main.c kseq.h hopo_counter.c hopo_counter.h hopo_set.c hopo_set.h
tatajuba_LDADD = $(LDADD)
