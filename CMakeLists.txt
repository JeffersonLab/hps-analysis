cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
cmake_policy(SET CMP0005 NEW)
cmake_policy(SET CMP0042 NEW)

project(hps-analysis)

add_subdirectory(MiniDst)
add_subdirectory(BumpHunt2016)
add_subdirectory(Moller)
add_subdirectory(EventMixer)
