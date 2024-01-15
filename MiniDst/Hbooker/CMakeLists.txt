# Create the MiniDst library
project(Hbooker)

ROOT_GENERATE_DICTIONARY(HbookerDict HBooker.h MODULE ${PROJECT_NAME} LINKDEF HBookerLinkDef.h )

add_library(${PROJECT_NAME} SHARED HBooker.cxx HBookerDict.cxx)
target_include_directories(${PROJECT_NAME} BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "HBooker.h")
# Selective linking:
#target_link_libraries(Dst2016 PUBLIC BaseAna ROOT::Core ROOT::Hist ROOT::Tree ROOT::Physics )
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR} )

# install ROOT pcm file
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
