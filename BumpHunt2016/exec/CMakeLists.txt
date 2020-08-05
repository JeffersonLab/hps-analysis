project(exec)

add_executable(run_bumphunt runbumphunt.cpp)
target_link_libraries(run_bumphunt BumpHunt2016 ${ROOT_LIBRARIES})
install(TARGETS run_bumphunt)