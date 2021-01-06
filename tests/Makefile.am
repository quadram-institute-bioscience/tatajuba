AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/src  -I$(top_srcdir)/submodules/bwa -I$(top_srcdir)/@BIOMCMCLIB@/lib \
							@OPENMP_CPPFLAGS@  @ZLIB_CFLAGS@ @LZMA_CFLAGS@ -DTEST_FILE_DIR='"$(srcdir)/files/"'
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(top_builddir)/submodules/bwa/libbwa.la $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la \
				$(GTKDEPS_LIBS) @CHECK_LIBS@  @ZLIB_LIBS@ @LZMA_LIBS@ $(AM_LDFLAGS)

EXTRA_DIST = files # directory with fasta etc files (accessed with #define TEST_FILE_DIR above)
# list of programs to be compiled only with 'make check':
check_PROGRAMS = check_bwa debug_bwa debug_bwamem
# list of test programs (duplicate, since we want all to be compiled only with 'make check'):
TESTS = $(check_PROGRAMS)           
