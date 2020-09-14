AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  -I$(top_srcdir)/@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@ @ZLIB_LIBS@
AM_CFLAGS = @SIMD_FLAGS@ @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) @ZLIB_LIBS@ $(top_builddir)/@BIOMCMCLIB@/lib/libbiomcmc_static.la  $(AM_LDFLAGS)

common_headers = kalign.h \
tldevel.h \
rng.h \
global.h \
alignment_parameters.h \
bisectingKmeans.h \
sequence_distance.h \
alignment.h \
bpm.h

common_src = run_kalign.c \
tldevel.c \
rng.c \
alignment_parameters.c \
bisectingKmeans.c \
sequence_distance.c \
alignment.c \
bpm.c

noinst_LTLIBRARIES = libkalign.la   ## change to noinst_LT... if should be linked statically (not installed globally)
libkalign_la_SOURCES = config.h $(common_headers) $(common_src) 

