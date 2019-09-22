cmake_minimum_required(VERSION 3.12)
project(oicugl VERSION 0.1.0)

#
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
endif()

#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/bin)

#
add_subdirectory(src/)
add_subdirectory(samples/)

# EOF