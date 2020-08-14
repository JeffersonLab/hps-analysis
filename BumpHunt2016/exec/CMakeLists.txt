project(exec)

add_executable(run_bumphunt runbumphunt.cpp)
target_link_libraries(run_bumphunt BumpHunt2016 ${ROOT_LIBRARIES})

#
# Somehow on Linux the stdc++ does not (always) get properly picked up.
#
IF(CMAKE_SYSTEM_NAME STREQUAL Linux)
   target_link_libraries(run_bumphunt -lstdc++ )   
endif(CMAKE_SYSTEM_NAME STREQUAL Linux)



install(TARGETS run_bumphunt)
