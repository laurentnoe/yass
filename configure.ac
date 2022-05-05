dnl Process this file with autoconf to produce a configure script
AC_PREREQ([2.60])
AC_INIT([yass],[1.16],[laurent.noe@univ-lille.fr])
AC_CONFIG_AUX_DIR(build)
AC_CANONICAL_TARGET
AC_PROG_MAKE_SET
AC_CONFIG_SRCDIR([src/main.c])

dnl Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([stdlib.h string.h assert.h math.h errno.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

dnl Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset clock pow sqrt strchr strdup strtol strtoul])


dnl 1) threads options

dnl abc) : with-threads option
AC_ARG_WITH(
	threads,
	[  --with-threads          compile with all threads],
	[threads="yes"],
	[threads="no"])

if test "$threads" = "yes"; then
  AC_MSG_RESULT(detected cpu parameter: threads);
 
  dnl pthread library here ? (UNIX system)
  AC_MSG_CHECKING(for pthread lib)
  AC_CHECK_LIB(pthread, pthread_create,
	[have_pthread="yes"],
	[have_pthread="no"])
  if test "$have_pthread" = "yes"; then
    AC_MSG_RESULT(yes)
    AC_DEFINE([HAVE_PTHREAD],[1],[pthread available on the computer])
       CFLAGS=" $CFLAGS -DTHREAD_ASSEMBLE_ALIGN -DTHREAD_FORWARD_REVERSE -DTHREAD_QUERY_CHUNK "
       LIBS="$LIBS -lpthread"
  else
    AC_MSG_RESULT(no : disabling \"--with-threads\" parameter)
  fi

else
   if test -s /proc/cpuinfo; then
      AC_MSG_CHECKING(for multi-processor)
      dnl Multithread advice
      if test `grep -c '^processor' /proc/cpuinfo` -gt 2; then
         AC_MSG_RESULT(yes) 
         AC_MSG_RESULT(- try \"configure --with-threads\")
      else
         AC_MSG_RESULT(no)
      fi;
   fi;
fi





dnl a) : with-thread-fr option
AC_ARG_WITH(
	thread-fr,
	[  --with-thread-fr        compile with two separate threads for Forward and Reverse sequence],
	[thread_forward_reverse="yes"],
	[thread_forward_reverse="no"])

if test "$thread_forward_reverse" = "yes"; then
  AC_MSG_RESULT(detected cpu parameter: thread-fr);
 
  dnl pthread library here ? (UNIX system)
  AC_MSG_CHECKING(for pthread lib)
  AC_CHECK_LIB(pthread, pthread_create,
	[have_pthread="yes"],
	[have_pthread="no"])
  if test "$have_pthread" = "yes"; then
    AC_MSG_RESULT(yes)
    AC_DEFINE([HAVE_PTHREAD],[1],[pthread available on the computer])
       CFLAGS=" $CFLAGS -DTHREAD_FORWARD_REVERSE "
       LIBS="$LIBS -lpthread"  
  else
    AC_MSG_RESULT(no : disabling \"--with-thread-fr\" parameter)
  fi

else
   if test -s /proc/cpuinfo; then
      AC_MSG_CHECKING(for multi-processor)
      dnl Multithread advice
      if test `grep -c '^processor' /proc/cpuinfo` -gt 1; then
         AC_MSG_RESULT(yes) 
         AC_MSG_RESULT(- try \"configure --with-thread-fr\")
      else
         AC_MSG_RESULT(no)
      fi;
   fi;
fi





dnl b) : with-thread-aa option
AC_ARG_WITH(
	thread-aa,
	[  --with-thread-aa        compile with two separate threads for Assemble and Align steps],
	[thread_assemble_align="yes"],
	[thread_assemble_align="no"])

if test "$thread_assemble_align" = "yes"; then
  AC_MSG_RESULT(detected cpu parameter: thread-aa);
 
  dnl pthread library here ? (UNIX system)
  AC_MSG_CHECKING(for pthread lib)
  AC_CHECK_LIB(pthread, pthread_create,
	[have_pthread="yes"],
	[have_pthread="no"])
  if test "$have_pthread" = "yes"; then
    AC_MSG_RESULT(yes)
    AC_DEFINE([HAVE_PTHREAD],[1],[pthread available on the computer])
       CFLAGS=" $CFLAGS -DTHREAD_ASSEMBLE_ALIGN "
       LIBS="$LIBS -lpthread"  
  else
    AC_MSG_RESULT(no : disabling \"--with-thread-aa\" parameter)
  fi

else
   if test -s /proc/cpuinfo; then
      AC_MSG_CHECKING(for multi-processor)
      dnl Multithread advice
      if test `grep -c '^processor' /proc/cpuinfo` -gt 1; then
         AC_MSG_RESULT(yes) 
         AC_MSG_RESULT(- try \"configure --with-thread-aa\")
      else
         AC_MSG_RESULT(no)
      fi;
   fi;
fi


dnl c) : with-thread-qc option
AC_ARG_WITH(
	thread-qc,
	[  --with-thread-qc        compile with threads for the query chunks],
	[thread_query_chunk="yes"],
	[thread_query_chunk="no"])

if test "$thread_query_chunk" = "yes"; then
  AC_MSG_RESULT(detected cpu parameter: thread-qc);
 
  dnl pthread library here ? (UNIX system)
  AC_MSG_CHECKING(for pthread lib)
  AC_CHECK_LIB(pthread, pthread_create,
	[have_pthread="yes"],
	[have_pthread="no"])
  if test "$have_pthread" = "yes"; then
    AC_MSG_RESULT(yes)
    AC_DEFINE([HAVE_PTHREAD],[1],[pthread available on the computer])
       CFLAGS=" $CFLAGS -DTHREAD_QUERY_CHUNK "
       LIBS="$LIBS -lpthread"  
  else
    AC_MSG_RESULT(no : disabling \"--with-thread-qc\" parameter)
  fi

else
   if test -s /proc/cpuinfo; then
      AC_MSG_CHECKING(for multi-processor)
      dnl Multithread advice
      if test `grep -c '^processor' /proc/cpuinfo` -gt 1; then
         AC_MSG_RESULT(yes) 
         AC_MSG_RESULT(- try \"configure --with-thread-qc\")
      else
         AC_MSG_RESULT(no)
      fi;
   fi;
fi


CFLAGS="$CFLAGS -O3 -Wall -funroll-loops -pipe -fomit-frame-pointer "
LIBS="$LIBS -lm"

AC_SUBST(CFLAGS)
AC_SUBST(LIBS)
AM_WITH_DMALLOC
AM_INIT_AUTOMAKE
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 src/Makefile])
AC_OUTPUT





