# Create the CERNLIB HBook maker executable.
# Because the PAWC common block must be in the lower 32-bits of memory space, this does not seem to work
# as a library, at least not with the 2005 version of CERNLIB which is currently used at JLab.
# It seems that only with an executable the memory can be put in the desired location.
# The error when trying this as a library is:
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# LOCB/LOCF: address 0x7f754d4abf60 exceeds the 32 bit address space
# or is not in the data segments
# This may result in program crash or incorrect results
# Therefore we will stop here
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# I thus only make the executable.
#
# A second issue is that the version of libgfortran.a needs to match the one used to
# compile the CERNLIB libraries. This library seems to be /usr/lib64/libgfortran.so.3.0.0
# This is setup in the FindCERLIB.cmake script in MiniDST/cmake
#
project(Hbooker)

# For the include of the MiniDst.h header file.
get_property(MiniDstLib_Include TARGET MiniDst PROPERTY INCLUDE_DIRECTORIES)
include_directories(BEFORE ${MiniDstLib_Include})

# No need for the dictionary if we are not making a library.
#
# ROOT_GENERATE_DICTIONARY(HbookerDict Hbooker.h MODULE ${PROJECT_NAME} LINKDEF HbookerLinkDef.h )
# add_library(${PROJECT_NAME} SHARED Hbooker.cxx HbookerDict.cxx)
#
add_executable(make_hbook_dst make_hbook_dst.cxx Hbooker.cxx)
target_include_directories(make_hbook_dst BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(make_hbook_dst PUBLIC MiniDst ${ROOT_LIBRARIES} ${CERNLIB_LIBRARIES})
if(APPLE)
target_link_options(make_hbook_dst PUBLIC "-Wl,-no_compact_unwind")
endif()
target_compile_options(make_hbook_dst PUBLIC "-Dg77Fortran" )
target_compile_options(make_hbook_dst PUBLIC "-no-pie")

# No need for a public header.
# set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "Hbooker.h")
# install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR} )
