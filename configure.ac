AC_INIT([tatajuba], [1.0.0], [leomrtns@gmail.com],[], [https://github.com/quadram-institute-bioscience/tatajuba])
AC_CONFIG_HEADER([biomcmc-lib/lib/config.h])
AC_CONFIG_AUX_DIR([.])

AC_CONFIG_SRCDIR([src/main.c])
AC_PREFIX_DEFAULT(${HOME})
dnl libtoolize recomendation
AC_CONFIG_MACRO_DIR([m4])
CFLAGS=${CFLAGS-""}  # Override default O2
CXXFLAGS=${CXXFLAGS-""}

AM_INIT_AUTOMAKE([-Wall -Werror]) # automake initialization (completely unrelated to compiler arguments)
AC_USE_SYSTEM_EXTENSIONS
AM_MAINTAINER_MODE
dnl must be called after AC_USE_SYSTEM_EXTENSIONS (new restriction in autoconf2.68) 
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_LANG_C
AC_PROG_CC
AC_PROG_INSTALL
AC_GNU_SOURCE # directs the libc header files to provide the standard GNU system interfaces including all GNU extensions
AC_PROG_LIBTOOL # Libtool (for library): in lib/Makefile.am, we include an "abstract" libfoo.la

# M4 macros for checking of CPU features from kalign3
m4_include([m4/ax_gcc_x86_avx_xgetbv.m4])
m4_include([m4/ax_gcc_x86_cpuid.m4])
m4_include([m4/ax_check_compile_flag.m4])
m4_include([m4/ax_ext.m4])
m4_include([m4/ax_openmp.m4])
AX_EXT

dnl openMP: sets $OPENMP_CFLAGS which should be passed to CFLAGS, CPPFLAGS; creates preprocessor macro _OPENMP 
dnl         (checked with "ifdef _OPENMP"); user can disable it through "--disable-openmp"
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CPPFLAGS)

PKG_CHECK_MODULES([CHECK], [check >= 0.9.10],  # unit tests in for "make check"
                  [AC_MSG_NOTICE([testing library `check` found])],
                  [AC_MSG_NOTICE([library `check` NOT found, you won't be able to run the battery of tests])]
                  )
PKG_CHECK_MODULES([ZLIB], [zlib],
                  AC_DEFINE([HAVE_ZLIB],[1],[]),
                  AC_DEFINE([HAVE_ZLIB],[0],[set to one if gz compression library found])
                  )
PKG_CHECK_MODULES([LZMA], [liblzma],
                  AC_DEFINE([HAVE_LZMA],[1],[]),
                  AC_DEFINE([HAVE_LZMA],[0],[set to one if xz compression library found])
                  )
AC_CHECK_HEADER(bzlib.h,
                  AC_DEFINE([HAVE_BZIP2],[1],[]),
                  AC_DEFINE([HAVE_BZIP2],[0],[set to one if bz2 compression library found])
                  )
AC_SEARCH_LIBS([BZ2_bzlibVersion], [bz2]) #  unlike AC_CHECK_LIB, does not define HAVE_LIBBZ2 (thus our custom define above)
AC_SEARCH_LIBS([sqrt], [m])

AC_MSG_RESULT([ ===  configuration options specific to $PACKAGE_NAME]) # library checks (math etc.) are done by biomcmc-lib

AC_MSG_CHECKING([whether to build (slower) debug code])
AC_ARG_ENABLE([debug], [AS_HELP_STRING([--enable-debug], [enable debugging with gdb and friends (default=no)])],[debugit="$enableval"], [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AC_DEFINE([BIOMCMC_DEBUG],[],[Debug Mode, with assert()-like error checking])
    AM_CFLAGS="${AM_CFLAGS} -g -W -pg -Wall -Werror -Wno-uninitialized -O0 -std=gnu11"
else
    AM_CFLAGS="${AM_CFLAGS} -funroll-loops -fomit-frame-pointer -finline-functions -O4 -std=gnu11"
fi

AC_MSG_CHECKING([whether you want static binaries])
AC_ARG_ENABLE(static-binary,
    [  --enable-static-binary  static binaries, that run on same arch without the libraries [[default=no]]],
    [ statbin_use="yes" ], [ statbin_use=""    ])
if test -n "${statbin_use}"; then
    AC_MSG_RESULT([yes])
    AM_LDFLAGS="-static ${AM_LDFLAGS}";
else
    AC_MSG_RESULT([no])
fi

AC_MSG_RESULT([ ===  end of specific configuration options])

dnl propagate changed vars among final makefiles
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([MPI_CXXLIBS])
AC_SUBST([MPI_CXXFLAGS])

dnl better solution is to allow arbitrary location for biomcmc-lib (it relies on ln)
AC_CHECK_FILE([${srcdir}/biomcmc-lib],[],[ln -s submodules/biomcmc-lib ${srcdir}/biomcmc-lib])
AC_CHECK_FILE([${srcdir}/biomcmc-lib/configure.ac],[], [AC_MSG_ERROR(["biomcmc-lib submodule missing, please git clone --recursive or link by hand to location of source code"])])
AC_SUBST([BIOMCMCLIB], [biomcmc-lib])
dnl generate makefiles (last recipes, after defining CFLAGS etc.)
AC_CONFIG_FILES([biomcmc-lib/Makefile biomcmc-lib/lib/Makefile biomcmc-lib/tests/Makefile 
                Makefile submodules/bwa/Makefile src/Makefile kalign/Makefile tests/Makefile])
AC_OUTPUT

echo \
"---- 
Configuration parameters for $PACKAGE_NAME-$PACKAGE_VERSION:

 Source code location:   ${srcdir}
 Compiler     (library): ${CC}
 Compiler flags: ${AM_CFLAGS} 
 Linker  flags:  ${AM_LDFLAGS}
 Install path:   ${prefix}
----"
