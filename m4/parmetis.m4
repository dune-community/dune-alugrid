# searches for PARMETIS header and lib 

AC_DEFUN([ALUGRID_PATH_PARMETIS],[
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([ALUGRID_PATH_METIS])
  AC_REQUIRE([DUNE_MPI])

  AC_ARG_WITH(parmetis,
    AC_HELP_STRING([--with-parmetis=PATH],[directory with ParMETIS (version >= 3.1)  inside]))

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"
LIBS=""

## do nothing if no --with-parmetis was supplied
if test x$with_parmetis != x && test x$with_parmetis != xno ; then

  if test "x$with_metis" = "xno" ; then
    AC_MSG_WARN([ParMETIS only works with METIS!])
  fi

  if test x$with_parmetis == xyes ; then
    AC_MSG_ERROR([You have to provide a directory --with-parmetis=PATH])
  fi

  if test -d $with_parmetis; then
    # expand tilde / other stuff
    PARMETISROOT=`cd $with_parmetis && pwd`
  else
    AC_MSG_ERROR([Path $with_parmetis supplied for --with-parmetis does not exist!])
  fi

PARMETIS_LIB_PATH="$PARMETISROOT/lib"
PARMETIS_INCLUDE_PATH="$PARMETISROOT/include"

if ! test -d $PARMETIS_LIB_PATH ; then 
  PARMETIS_LIB_PATH="$PARMETISROOT"
fi 
if ! test -d $PARMETIS_INCLUDE_PATH ; then 
  PARMETIS_INCLUDE_PATH="$PARMETISROOT"
fi 

# set variables so that tests can use them
REM_CPPFLAGS=$CPPFLAGS

LDFLAGS="$LDFLAGS -L$PARMETIS_LIB_PATH -lm"
CPPFLAGS="$CPPFLAGS $PARMETISDEF -I$PARMETIS_INCLUDE_PATH"

# check for header
  AC_LANG_PUSH([C])
  REM_CC="$CC"
  CC="$MPICC"
  AC_CHECK_HEADER([parmetis.h], 
   [PARMETIS_CPPFLAGS="$PARMETISDEF -I$PARMETIS_INCLUDE_PATH"
	HAVE_PARMETIS="1"],
  AC_MSG_WARN([parmetis.h not found in $PARMETIS_INCLUDE_PATH]))

CPPFLAGS="$REM_CPPFLAGS $PARMETISDEF -I$PARMETIS_INCLUDE_PATH"
REM_CPPFLAGS=

REM_LDFLAGS=$LDFLAGS
# add math library 
LDFLAGS="$LDFLAGS"
LIBS="-lmetis -lm"

# if header is found...
if test x$HAVE_PARMETIS = x1 ; then
  AC_CHECK_LIB(parmetis,[ParMETIS_V3_PartKway],
	[PARMETIS_LIBS="-lparmetis"
         PARMETIS_LDFLAGS="-L$PARMETIS_LIB_PATH"
         LIBS="$LIBS $PARMETIS_LIBS"],
	[HAVE_PARMETIS="0"
	AC_MSG_WARN(libparmetis.a not found in $PARMETIS_LIB_PATH)])
fi

LDFLAGS=$REM_LDFLAGS

# pop default language 
CC="$REM_CC"
AC_LANG_POP([C])

## end of parmetis check (--without wasn't set)
fi

# survived all tests?
if test x$HAVE_PARMETIS = x1 ; then
  AC_SUBST(PARMETIS_LIBS, $PARMETIS_LIBS)
  AC_SUBST(PARMETIS_LDFLAGS, $PARMETIS_LDFLAGS)
  AC_SUBST(PARMETIS_CPPFLAGS, $PARMETIS_CPPFLAGS)
  AC_DEFINE(HAVE_PARMETIS, 1, [Define to 1 if parmetis-library is found])

  
  # add to global list
  DUNE_ADD_ALL_PKG([PARMETIS], [\${PARMETIS_CPPFLAGS}],
                   [\${PARMETIS_LDFLAGS}], [\${PARMETIS_LIBS}])

  ALU_PARMETIS_VFILE=$PARMETISROOT/VERSION
  ALU_PARMETIS_VERSION="(unknown)"
  if test -f $PARMETISROOT/VERSION ; then 
    ALU_PARMETIS_VERSION="(Version `cat $ALU_PARMETIS_VFILE | head -1 | cut -d " " -f 3`)"
  elif test -f $PARMETIS_INCLUDE_PATH/parmetis.h ; then 
    PARMETIS_HEADER=$PARMETIS_INCLUDE_PATH/parmetis.h
    VERSION_MAJOR=`grep "PARMETIS_MAJOR_VERSION" $PARMETIS_HEADER | awk '{ print [$]3 }'`
    VERSION_MINOR=`grep "PARMETIS_MINOR_VERSION" $PARMETIS_HEADER | awk '{ print [$]3 }'`
    VERSION_SUBMINOR=`grep "PARMETIS_SUBMINOR_VERSION" $PARMETIS_HEADER | awk '{ print [$]3 }'`
    ALU_PARMETIS_VERSION="($VERSION_MAJOR.$VERSION_MINOR.$VERSION_SUBMINOR)"
  fi  

  # set variable for summary
  with_parmetis="yes $ALU_PARMETIS_VERSION"
else
  AC_SUBST(PARMETIS_LIBS, "")
  AC_SUBST(PARMETIS_LDFLAGS, "")
  AC_SUBST(PARMETIS_CPPFLAGS, "")

  # set variable for summary
  with_parmetis="no"
fi
  
# also tell automake
AM_CONDITIONAL(PARMETIS, test x$HAVE_PARMETIS = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

])
