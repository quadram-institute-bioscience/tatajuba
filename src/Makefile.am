AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/submodules/bwa -I$(top_srcdir)/@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@  
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@
LDADD = $(top_builddir)/submodules/bwa/libbwa.la $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la \
				$(GTKDEPS_LIBS) @ZLIB_LIBS@ @LZMA_LIBS@ $(AM_LDFLAGS)

bin_PROGRAMS = tatajuba 
tatajuba_SOURCES = main.c kseq.h genome_set.c genome_set.h context_histogram.c context_histogram.h hopo_counter.c hopo_counter.h
tatajuba_LDADD = $(LDADD)
