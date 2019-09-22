cmake_minimum_required(VERSION 3.12)
project(oicugl VERSION 0.1.0)

#
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
endif()

#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

#
add_subdirectory(src/)

# EOF