AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../src -I$(srcdir)/../kalign -I$(srcdir)/../bwa \
							-I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ -DTEST_FILE_DIR='"$(srcdir)/files/"'
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  ../kalign/libkalign.la ../bwa/libbwa.la ../biomcmc-lib/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

EXTRA_DIST = files # directory with fasta etc files (accessed with #define TEST_FILE_DIR above)
# list of programs to be compiled only with 'make check':
check_PROGRAMS = check_bwa debug_bwa
# list of test programs (duplicate, since we want all to be compiled only with 'make check'):
TESTS = $(check_PROGRAMS)           
