AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/kalign -I$(top_srcdir)/submodules/bwa -I$(top_srcdir)/@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  $(top_builddir)/kalign/libkalign.la $(top_builddir)/submodules/bwa/libbwa.la $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

bin_PROGRAMS = tatajuba 
tatajuba_SOURCES = main.c kseq.h genome_set.c genome_set.h  context_histogram.c context_histogram.h
tatajuba_LDADD = $(LDADD)

BUILT_SOURCES =  histogram_comparison.c
