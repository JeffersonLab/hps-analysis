# Just in case I don't gripe about Hpstr enough. Naming the library "event" is a bad idea.
# It is a common name and will conflict with other libraries. For instance, the event library for system events.
# But, that is what we got, so now we have to force cmake to find *the correct* event library.

SET(HPSTR_DIR "/" CACHE STRING "Directory pointing to the HPSTR sources.")
message(STATUS "HPSTR_DIR: ${HPSTR_DIR}")
find_path(HPSTR_INCLUDE_DIR CalCluster.h HINT  /data/HPS/hpst/event/include ${HPSTR_DIR}/event/include)
#set(HPSTR_INCLUDE_DIRS ${HPSTR_INCLUDE_DIR})
message(STATUS "HPSTR_INCLUDE_DIR: ${HPSTR_INCLUDE_DIR}")

function(event_library_check validator_result_var item)

    # message(STATUS "==========++++++++++++++++++++++++ Checking for ${item} ++++++++++++++++++++++++++")

    file(WRITE ${CMAKE_BINARY_DIR}/event_library_check.cxx
        "#include <CalCluster.h>
         int main(){
         CalCluster *clust = new CalCluster();
          }\n"
    )
    #message(STATUS "HPSTR_INCLUDE_DIR: ${HPSTR_INCLUDE_DIR}")
    set(CMAKE_CXX_FLAGS "-I${HPSTR_INCLUDE_DIR} -I${ROOT_INCLUDE_DIRS} ${CMAKE_CXX_FLAGS}")
    #message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

    # -DINCLUDE_DIRECTORIES=/data/HPS/hpstr/event/include
    try_compile(HAS_SYMBOL SOURCES ${CMAKE_BINARY_DIR}/event_library_check.cxx
            CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}  LINK_LIBRARIES ${item} ${ROOT_LIBRARIES} OUTPUT_VARIABLE OUTPUT)
    if(HAS_SYMBOL)
        #message(STATUS "We seem to be able to compile and link against ${item}")
        #message(STATUS "Output: ${OUTPUT}")
        set(${validator_result_var} TRUE PARENT_SCOPE)
    else()
        #message(STATUS "We seem to NOT be able to compile and link against ${item}")
        #message(STATUS "Output: ${OUTPUT}")
        set(${validator_result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()


find_library(HPSTR_LIBRARY event HINT ${HPSTR_LIBDIR} VALIDATOR event_library_check)
#find_library(HPSTR_LIBRARY event HINT ${HPSTR_LIBDIR} )
#set(HPSTR_LIBRARIES ${HPSTR_LIBRARY})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HPSTR DEFAULT_MSG HPSTR_LIBRARY HPSTR_INCLUDE_DIR)

if( HPSTR_FOUND)
    set(HPSTR_LIBRARIES ${HPSTR_LIBRARY})
    set(HPSTR_INCLUDE_DIRS ${HPSTR_INCLUDE_DIR})
    message(STATUS "HPSTR found: ${HPSTR_INCLUDE_DIRS} ${HPSTR_LIBRARIES}")
else()
    message(STATUS "==================================================================")
    message(STATUS "Hpstr was not found, so you will not be able to convert Hpstr files.")
    message(STATUS "Sorry. ")
    message(STATUS "You really should not have to convert Hpstr files, if that format was sane. ")
    message(STATUS "It is not. Go and ask the creators and promotors of that format to explain it to you.")
    message(STATUS "While you are talking to them, also ask that they fix their CMakeLists.txt to properly install stuff.")
    message(STATUS "In the meantime, just convert the SLCIO files directly to minidst.")
    message(STATUS "=================================================================")
endif(HPSTR_FOUND)