cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0042 NEW)

project(hps-analysis)

add_subdirectory(MiniDst)
add_subdirectory(BumpHunt2016)
add_subdirectory(Moller)
add_subdirectory(EventMixer)
add_subdirectory(Utilities)

# Add the uninstall target
if(NOT TARGET uninstall)
    configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
            IMMEDIATE @ONLY)

    add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
