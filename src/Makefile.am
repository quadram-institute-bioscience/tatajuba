AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../kalign -I$(srcdir)/../bwa -I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  ../kalign/libkalign.la ../bwa/libbwa.la ../biomcmc-lib/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = tatajuba 
tatajuba_SOURCES = main.c kseq.h context_histogram.c genome_set.c context_histogram.h genome_set.h 
tatajuba_LDADD = $(LDADD)
