# Create the Dst2016 library
cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
project(Dst2016)

# link_directories("/data/HPS/lib")

# The ROOT_GENERATE_DICTIONARY does not seem to respect the included target_link_libraries target.
# We thus need to get the include directories from the target by hand and add them to the include path.
# These need to go at the END of this file, so that the BEFORE statememt properly puts the include directories early.
# These directories should be listed before the location where the headers are installed, so that changes to header files
# are properly picked up.
get_property(BaseAna_Include TARGET BaseAna2 PROPERTY INCLUDE_DIRECTORIES)
include_directories(BEFORE ${BaseAna_Include})

get_property(MiniDstLib_Include TARGET MiniDst PROPERTY INCLUDE_DIRECTORIES)
include_directories(BEFORE ${MiniDstLib_Include})

ROOT_GENERATE_DICTIONARY(Dst2016Dict Dst2016.h MODULE ${PROJECT_NAME} LINKDEF Dst2016LinkDef.h )

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
