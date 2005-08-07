# check for thing needed to compile the ALU-Grid library. 

AC_DEFUN([ALU3D_CHECK_ALL],[
  AC_LANG_PUSH([C++])
dnl check for programs
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([AC_PROG_CXXCPP])

  AC_REQUIRE([AC_PROG_INSTALL])
  AC_REQUIRE([AC_PROG_LN_S])
  AC_REQUIRE([AC_PROG_MAKE_SET])
  AC_REQUIRE([AC_PROG_RANLIB])
  AC_REQUIRE([AC_PROG_LIBTOOL])

dnl checks for header files.
  AC_REQUIRE([AC_HEADER_STDC])
  AC_CHECK_HEADERS([malloc.h string.h])

dnl checks for typedefs, structures, and compiler characteristics.
  AC_REQUIRE([AC_C_CONST])
  AC_REQUIRE([AC_C_INLINE])
  AC_REQUIRE([AC_TYPE_SIZE_T])
  AC_REQUIRE([AC_STRUCT_TM])

dnl check for library functions
  AC_REQUIRE([AC_FUNC_MALLOC])
  AC_LANG_PUSH([C++])
  AC_CHECK_LIB([m], [pow])
  AC_CHECK_FUNCS([sqrt strchr])
  AC_LANG_POP

dnl check all components
  AC_REQUIRE([ALU3D_PATH_METIS])
  AC_REQUIRE([ALU3D_PATH_PARTY])
  AC_REQUIRE([ALU3D_SERIAL_PARALLEL])

  # convenience-variables if every found package should be used
  AC_SUBST(ALL_PKG_LIBS, "$LIBS $ALU3D_PKG_LIBS")
  AC_SUBST(ALL_PKG_LDFLAGS, "$LDFLAGS $ALU3D_PKG_LDFLAGS")
  AC_SUBST(ALL_PKG_CPPFLAGS, "$CPPFLAGS $ALU3D_PKG_CPPFLAGS")
])


AC_DEFUN([ALU3D_SUMMARY_ALL],[
  # show search results
  echo
if test x$with_parallel = xnone ; then
  echo "ALU3d-Grid is not ready to compile for computations!"
else
  echo "ALU3d-Grid is ready to compile for $with_parallel computations!"
fi

if test x$with_parallel = xparallel ; then
  echo
  echo The following components where found: 
  echo "------------------------------------"
  echo  
  echo "METIS............: $with_metis"
  echo "PARTY............: $with_party"
  echo "MPI..............: $with_mpi"
  echo
  echo "------------------------------------"
fi
  echo
  echo "See ./configure --help and config.log for reasons why a component wasn't found"
  echo

])
