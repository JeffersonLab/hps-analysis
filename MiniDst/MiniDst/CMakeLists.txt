# Create the MiniDst library
project(MiniDst)

ROOT_GENERATE_DICTIONARY(MiniDstDict MiniDst.h MODULE ${PROJECT_NAME} LINKDEF MiniDstLinkDef.h )

add_library(${PROJECT_NAME} SHARED MiniDst.cxx MiniDstDict.cxx)
target_include_directories(${PROJECT_NAME} BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "MiniDst.h")
# Selective linking:
#target_link_libraries(Dst2016 PUBLIC BaseAna ROOT::Core ROOT::Hist ROOT::Tree ROOT::Physics )
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR} )

# install ROOT pcm file
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
