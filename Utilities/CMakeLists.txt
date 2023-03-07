cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0042 NEW)

project(Utilties VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

#
# Seem you need to set these (at least both library ones) to get the rdict.pcm and .rootmap files
# to be installed properly. (?)
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin" CACHE PATH "Build location of binary executatbles.")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" CACHE PATH "Build location of shared libraries.")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib" CACHE PATH "Build location of static libraries.")

set(CMAKE_INSTALL_LIBDIR           "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation location of libraries. ")
set(CMAKE_INSTALL_BINDIR           "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation location of binaries. ")
set(CMAKE_INSTALL_INCDIR           "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Installation location of header files. ")
set(CMAKE_INSTALL_DOCDIR           "${CMAKE_INSTALL_PREFIX}/doc" CACHE PATH "Installation location of documentation. ")

# list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame ROOTVecOps)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
message(STATUS "ROOT found at: ${ROOT_DIR} using use file: ${ROOT_USE_FILE}")

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

set(Utilities_PUBLIC_HEADERS "Utilities.h" "HPSMagField.h")
ROOT_GENERATE_DICTIONARY(UtilitiesDict ${Utilities_PUBLIC_HEADERS} MODULE ${PROJECT_NAME} LINKDEF UtilitiesLinkDef.h)

add_library(${PROJECT_NAME} SHARED Utilities.cxx HPSMagField.cxx UtilitiesDict.cxx)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${Utilities_PUBLIC_HEADERS}")

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR} )

link_directories("/data/HPS/lib")
