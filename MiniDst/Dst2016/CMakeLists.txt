# Create the Dst2016 library
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(Dst2016)

include_directories("../BaseAna")
#
# The ROOT_GENERATE_DICTIONARY does not seem to respect the included target_link_libraries target.
# We thus need to get the include directories from the target by hand and add them to the include path.
#
get_property(MiniDstLib_Include TARGET MiniDst PROPERTY INCLUDE_DIRECTORIES)
include_directories(${MiniDstLib_Include})

ROOT_GENERATE_DICTIONARY(Dst2016Dict Dst2016.h MODULE ${PROJECT_NAME} LINKDEF Dst2016LinkDef.h )

link_directories("/data/HPS/lib")

add_library(${PROJECT_NAME} SHARED Dst2016.cxx Dst2016Dict.cxx)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC MiniDst BaseAna2 ${ROOT_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "Dst2016.h")
# Selective linking:
#target_link_libraries(${PROJECT_NAME} PUBLIC BaseAna ROOT::Core ROOT::Hist ROOT::Tree ROOT::Physics )
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR})

# install ROOT pcm file
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
