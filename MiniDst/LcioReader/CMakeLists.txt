# Create the LcioReader library
project(LcioReader)

# The ROOT_GENERATE_DICTIONARY does not seem to respect the included target_link_libraries target.
# We thus need to get the include directories from the target by hand and add them to the include path.
# These need to go at the END of this file, so that the BEFORE statememt properly puts the include directories early.
# These directories should be listed before the location where the headers are installed, so that changes to header files
# are properly picked up.
get_property(MiniDstLib_Include TARGET MiniDst PROPERTY INCLUDE_DIRECTORIES)
include_directories(BEFORE ${MiniDstLib_Include})

ROOT_GENERATE_DICTIONARY(LcioReaderDict LcioReader.h MODULE ${PROJECT_NAME} LINKDEF LcioReaderLinkDef.h)

add_library(${PROJECT_NAME} SHARED LcioReader.cxx LcioReaderDict.cxx)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC MiniDst ${LCIO_LIBRARIES} ${ROOT_LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} )
