project(exec_moller)

add_executable(run_moller runmoller.cpp)
target_link_libraries(run_moller Moller ${ROOT_LIBRARIES})

#
# Somehow on Linux the stdc++ does not (always) get properly picked up.
#
IF(CMAKE_SYSTEM_NAME STREQUAL Linux)
   target_link_libraries(run_moller -lstdc++ )
endif(CMAKE_SYSTEM_NAME STREQUAL Linux)



install(TARGETS run_moller)
