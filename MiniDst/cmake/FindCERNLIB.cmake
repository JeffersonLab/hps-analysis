find_path(CERNLIB_CFORTRAN_INCLUDE_DIR cfortran.h ${CERN}/include/cfortran)
find_path(CERNLIB_INCLUDE_DIR hbook ${CERN}/include ${CERN}/include/cernlib/2023)

find_library(CERNLIB_LIBRARY_PACKLIB packlib ${CERN}/lib ${CERN}/lib/cernlib/2023/lib)
find_library(CERNLIB_LIBRARY_KERNLIB kernlib ${CERN}/lib ${CERN}/lib/cernlib/2023/lib)
set(CERNLIB_LIBRARIES ${CERNLIB_LIBRARY_PACKLIB} ${CERNLIB_LIBRARY_KERNLIB})
if(APPLE)
    set(CERNLIB_LIBRARIES ${CERNLIB_LIBRARIES} gfortran m dl)
else(APPLE)
    set(CERNLIB_LIBRARIES ${CERNLIB_LIBRARIES} gfortran m nsl crypt dl)
endif(APPLE)
#message("Variable CERN: ${CERN}")
#message("We found cernlib: ${CERNLIB_TEST}")
#message("CERNLIB_INCLUDE_DIR: ${CERNLIB_INCLUDE_DIR}")
#message("CERNLIB_LIBRARIES: ${CERNLIB_LIBRARIES}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(CERNLIB DEFAULT_MSG CERNLIB_LIBRARIES CERNLIB_INCLUDE_DIR )
