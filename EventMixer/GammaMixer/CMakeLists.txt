# Create the GammaMixer library
project(GammaMixer)

# The ROOT_GENERATE_DICTIONARY does not seem to respect the included target_link_libraries target.
# We thus need to get the include directories from the target by hand and add them to the include path.
# These need to go at the END of this file, so that the BEFORE statememt properly puts the include directories early.
# These directories should be listed before the location where the headers are installed, so that changes to header files
# are properly picked up.
get_property(MiniDstLib_Include TARGET MiniDst PROPERTY INCLUDE_DIRECTORIES)
include_directories(BEFORE ${MiniDstLib_Include})

ROOT_GENERATE_DICTIONARY(GammaMixerDict GammaMixer.h MODULE ${PROJECT_NAME} LINKDEF GammaMixerLinkDef.h )

add_library(${PROJECT_NAME} SHARED GammaMixer.cxx GammaMixerDict.cxx)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(${PROJECT_NAME} PUBLIC MiniDst ${ROOT_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "GammaMixer.h")
# Selective linking:
#target_link_libraries(Dst2016 PUBLIC BaseAna ROOT::Core ROOT::Hist ROOT::Tree ROOT::Physics )
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR} )

# install ROOT pcm file
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
