# Create the MiniDst library
project(Hbooker)

# For the include of the MiniDst.h header file.
get_property(MiniDstLib_Include TARGET MiniDst PROPERTY INCLUDE_DIRECTORIES)
include_directories(BEFORE ${MiniDstLib_Include})

ROOT_GENERATE_DICTIONARY(HbookerDict Hbooker.h MODULE ${PROJECT_NAME} LINKDEF HbookerLinkDef.h )

add_library(${PROJECT_NAME} SHARED Hbooker.cxx HbookerDict.cxx)
target_include_directories(${PROJECT_NAME} BEFORE PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(${PROJECT_NAME} PUBLIC MiniDst ${ROOT_LIBRARIES} ${CERNLIB_LIBRARIES})
if(APPLE)
target_link_options(${PROJECT_NAME} PUBLIC "-Wl,-no_compact_unwind")
endif()
target_compile_options(${PROJECT_NAME} PUBLIC "-Dg77Fortran" )
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "Hbooker.h")
# Selective linking:
#target_link_libraries(Dst2016 PUBLIC BaseAna ROOT::Core ROOT::Hist ROOT::Tree ROOT::Physics )
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCDIR} )

# install ROOT pcm file
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libDst2016.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
