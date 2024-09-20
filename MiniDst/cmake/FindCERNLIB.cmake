find_path(CERNLIB_CFORTRAN_INCLUDE_DIR cfortran.h ${CERN}/include/cfortran)
find_path(CERNLIB_INCLUDE_DIR hbook ${CERN}/include ${CERN}/include/cernlib/2023)

find_library(CERNLIB_LIBRARY_PACKLIB packlib ${CERN}/lib ${CERN}/lib/cernlib/2023/lib)
find_library(CERNLIB_LIBRARY_KERNLIB kernlib ${CERN}/lib ${CERN}/lib/cernlib/2023/lib)
set(CERNLIB_LIBRARIES ${CERNLIB_LIBRARY_PACKLIB} ${CERNLIB_LIBRARY_KERNLIB})
if(APPLE)
    # Sorry, this is not generic. It could be replaced with a proper "find" routine.
    set(CERNLIB_LIBRARIES ${CERNLIB_LIBRARIES} -L/opt/homebrew/Cellar/gcc/14.1.0_1/lib/gcc/current gfortran m dl)
else(APPLE)
    # We assume JLab. We need the *specific* version of the gfortran library.
    set(CERNLIB_LIBRARIES ${CERNLIB_LIBRARIES} /usr/lib64/libgfortran.so.3.0.0 m nsl crypt dl)
endif(APPLE)
#message("Variable CERN: ${CERN}")
#message("We found cernlib: ${CERNLIB_TEST}")
#message("CERNLIB_INCLUDE_DIR: ${CERNLIB_INCLUDE_DIR}")
#message("CERNLIB_LIBRARIES: ${CERNLIB_LIBRARIES}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(CERNLIB DEFAULT_MSG CERNLIB_LIBRARIES CERNLIB_INCLUDE_DIR )
