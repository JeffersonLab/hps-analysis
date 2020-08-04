cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(BaseAna2)

ROOT_GENERATE_DICTIONARY(BaseAnaDict BaseAna.h MODULE ${PROJECT_NAME} LINKDEF BaseAna_LinkDef.h )

add_library(${PROJECT_NAME} STATIC BaseAna.cpp BaseAnaDict.cxx
        HpsEvent.cxx EcalCluster.cxx SvtTrack.cxx GblTrack.cxx
        SvtHit.cxx HpsParticle.cxx EcalHit.cxx MCParticle.cxx HpsMCParticle.cxx )
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES})
set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

# Selective linking:
# target_link_libraries(${PROJECT_NAME} PUBLIC BaseAna ROOT::Core ROOT::Hist ROOT::Tree ROOT::Physics )
# install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR} )
