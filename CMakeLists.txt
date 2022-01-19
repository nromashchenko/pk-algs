cmake_minimum_required(VERSION 3.16)
project(xpas_algs)

if(EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    conan_basic_setup()
else()
    message(WARNING "The file conanbuildinfo.cmake doesn't exist, you have to run conan install first")
endif()

set(CMAKE_CXX_STANDARD 17)

add_executable(xpas_algs
        main.cpp
        common.cpp
        brute_force.cpp
        dc.cpp
        bb.cpp
        matrix.cpp
        ar.cpp)

target_link_libraries(xpas_algs ${CONAN_LIBS})


add_executable(test_matrix
        test_ranges.cpp
        matrix.cpp)
target_link_libraries(test_matrix ${CONAN_LIBS})

