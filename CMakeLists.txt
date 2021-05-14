cmake_minimum_required(VERSION 3.5)
project(oicugl VERSION 0.1.0)

#
set(CMAKE_CXX_STANDARD 14)

#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)
#
# Omit the next line to natively install CUGL via CMake
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/bin)

#
add_subdirectory(src/)
# Uncommment after a successful `make install`
#add_subdirectory(samples/)
# EOF
