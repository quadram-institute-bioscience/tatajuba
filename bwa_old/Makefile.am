AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(top_srcdir)/@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@  @ZLIB_LIBS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @CHECK_CFLAGS@ -lz -D_FILE_OFFSET_BITS=64 -fsigned-char -w 
LDADD = $(GTKDEPS_LIBS) @CHECK_LIBS@  $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la @ZLIB_LIBS@  $(AM_LDFLAGS)

original_headers = \
	bntseq.h \
	bwa.h \
	bwamem.h \
	bwase.h \
	bwtaln.h \
	bwtgap.h \
	bwt.h \
	bwt_lite.h \
	bwtsw2.h \
	kbtree.h \
	kseq.h \
	khash.h \
	ksort.h \
	ksort.h \
	kstring.h \
	ksw.h \
	kvec.h \
	bamlite.h \
	malloc_wrap.h \
	QSufSort.h \
	rle.h \
	rope.h \
	utils.h 
                 
original_src = \
	bwa.c \
	bamlite.c \
	bntseq.c \
	bwamem.c \
	bwamem_extra.c \
	bwamem_pair.c \
	bwape.c \
	bwase.c \
	bwaseqio.c \
	bwashm.c \
	bwtaln.c \
	bwt.c \
	bwtgap.c \
	bwt_gen.c \
	bwtindex.c \
	bwt_lite.c \
	bwtsw2_aux.c \
	bwtsw2_chain.c \
	bwtsw2_core.c \
	bwtsw2_main.c \
	bwtsw2_pair.c \
	fastmap.c \
	is.c \
	kopen.c \
	kthread.c \
	kstring.c \
	ksw.c \
	malloc_wrap.c \
	maxk.c \
	pemerge.c \
	QSufSort.c \
	rle.c \
	rope.c \
	utils.c 

noinst_LTLIBRARIES = libbwa.la # convenience lib (not installed globally)
libbwa_la_SOURCES = $(original_src) $(original_headers) wrapper_bwa.c wrapper_bwa.h 
