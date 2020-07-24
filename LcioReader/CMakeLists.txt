# Create the LcioReader library
project(LcioReader)

ROOT_GENERATE_DICTIONARY(LcioReaderDict LcioReader.h MODULE ${PROJECT_NAME} LINKDEF LcioReaderLinkDef.h)

add_library(${PROJECT_NAME} SHARED LcioReader.cxx LcioReaderDict.cxx)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} )
