#
# Specify the minimum CMake version required.  If the cmake version 
# minimum is not met, cmake will fail.
cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0042 NEW)

# Set the name of the project
project(MakeMiniDst)

# Set the project version
set(HPS_ANALYSIS_VERSION_MAJOR 0)
set(HPS_ANALYSIS_VERSION_MINOR 1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Offer the user the choice of overriding the installation directories
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-DDEBUG=1)
endif()

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) 
        SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "HPS_Analysis install prefix, default is project top dir." FORCE) 
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) 

message("The installation prefix is set to:  ${CMAKE_INSTALL_PREFIX} ")   

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin" CACHE PATH "Build location of binary executatbles.")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" CACHE PATH "Build location of shared libraries.")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" CACHE PATH "Build location of static libraries.")

set(CMAKE_INSTALL_LIBDIR           "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation location of libraries. ")
set(CMAKE_INSTALL_BINDIR           "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation location of binaries. ")
set(CMAKE_INSTALL_INCDIR           "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Installation location of header files. ")
set(CMAKE_INSTALL_DOCDIR           "${CMAKE_INSTALL_PREFIX}/doc" CACHE PATH "Installation location of documentation. ")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# You need to tell CMake where to find the ROOT installation.
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
# find ROOT installation
find_package(ROOT REQUIRED COMPONENTS Physics)
# COMPONENTS Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer
# Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame PyROOT Geom Gui
message(STATUS "ROOT found at: ${ROOT_DIR}")
message(STATUS " ")
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

find_package(LCIO)
if(LCIO_FOUND)
    include_directories( AFTER "${LCIO_INCLUDE_DIR}")
    message(STATUS "LCIO includes: ${LCIO_INCLUDE_DIR}")
else()
    message(STATUS "=================================================================")
    message(STATUS "LCIO subsystem was not found, so you will not be able to create  ")
    message(STATUS "root files from slcio file. To build LCIO go to MiniDST/lcio and ")
    message(STATUS "run: mkdir build; cd build; ")
    message(STATUS "cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} .. ")
    message(STATUS "=================================================================")
endif(LCIO_FOUND)

find_package(CERNLIB)
if(CERNLIB_FOUND)
    include_directories( AFTER "${CERNLIB_INCLUDE_DIR}")
    message(STATUS "CERNLIB found: ${CERNLIB_INCLUDE_DIR} ${CERNLIB_LIBRARIES}")
else()
    message(STATUS "==================================================================")
    message(STATUS "CERNLIB subsystem was not found, so you will not be able to create")
    message(STATUS "HBOOK output file. Sorry.")
    message(STATUS "To fix this set the variable CERN to you installation of cernlib. ")
    message(STATUS " - - - ")
    message(STATUS "=================================================================")
endif(CERNLIB_FOUND)


SET(CMAKE_INCLUDE_CURRENT_DIR ON)

if(APPLE)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

  SET(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    message("CMAKE_INSTALL_RPATH set to ${CMAKE_INSTALL_RPATH} ")
    
  ENDIF("${isSystemDir}" STREQUAL "-1")
endif(APPLE)

add_subdirectory(BaseAna)
add_subdirectory(MiniDst)
add_subdirectory(Dst2016)

if(LCIO_FOUND)
    add_subdirectory(LcioReader)
    add_executable(make_mini_dst make_mini_dst.cxx)
    target_link_libraries(make_mini_dst MiniDst Dst2016 LcioReader ${ROOT_LIBRARIES})
    install(TARGETS make_mini_dst)
    add_executable(lcio_test_reader lcio_test_reader.cpp)
    target_link_libraries(lcio_test_reader ${LCIO_LIBRARIES})
    install(TARGETS lcio_test_reader)

    if(CERNLIB_FOUND)
        add_subdirectory(Hbooker)
    endif(CERNLIB_FOUND)

endif(LCIO_FOUND)


# link_directories(${CMAKE_INSTALL_PREFIX}/lib)

#
# Somehow on Linux the stdc++ does not (always) get properly picked up.
#
#IF(CMAKE_SYSTEM_NAME STREQUAL Linux)
#   target_link_libraries(make_mini_dst -lstdc++ )
#   target_link_libraries(lcio_test_reader -lstdc++ )
#endif(CMAKE_SYSTEM_NAME STREQUAL Linux)

message("-----------------------------------------------")
message("* Configuration Summary: ")
message("*")
message("* ROOT version: ${ROOT_VERSION} ")
message("* ROOT location: $ENV{ROOTSYS} ")
message("* LCIO include: ${LCIO_INCLUDE_DIR} ")
message("* Compiler ID:  ${CMAKE_CXX_COMPILER_ID}")
message("* Compiler Ver: ${CMAKE_CXX_COMPILER_VERSION}")
