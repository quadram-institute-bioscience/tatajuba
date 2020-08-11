AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@ -lz -D_FILE_OFFSET_BITS=64 -fsigned-char -w 
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS)

common_headers =  bseq.h kalloc.h kdq.h ketopt.h khash.h kseq.h ksort.h ksw2.h kthread.h kvec.h minimap.h mmpriv.h sdust.h
                 
common_src =  align.c bseq.c chain.c esterr.c example.c format.c hit.c index.c kalloc.c ksw2_dispatch.c ksw2_extd2_sse.c \
		ksw2_exts2_sse.c ksw2_extz2_sse.c ksw2_ll_sse.c kthread.c main.c map.c misc.c options.c pe.c sdust.c sketch.c splitidx.c

noinst_LTLIBRARIES = libminimap2.la # convenience lib (not installed globally)
libminimap2_la_SOURCES = $(common_src) $(common_headers) 
