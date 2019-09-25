# FindCUGL.cmake


## Variables

1. [Include Directories] **CUGL_INCLUDE_DIR** 
2. [Link Directories] **CUGL_LIBRARY** 
3. [Boolean] **CUGL_INCLUDE_DIR_FOUND** 
4. [Boolean] **CUGL_LIBRARY_FOUND**
5. [Environment Variable Include] **CUGL_INCLUDE**
6. [Environment Variable Library] **CUGL_REDIST**

## Steps

#### Scientific Linux

1. Make sure you have CMake installed.
2. Clone the repository.
3. `cd CUGL`
4. `mkdir build && cd build`
5. `cmake .. -L`
6. `make`
7. `make install`
8. `setenv CUGL_INCLUDE PATH` where PATH is the absolute path to `$HOME/CUGL/build/bin/include`
9. `setenv CUGL_REDIST PATH` where PATH is the absolute path to `$HOME/CUGL/build/bin/lib`
10. Uncomment `add_subdirectory(samples)` in the root `CMakeLists.txt`
11. `make`
12. `.build/samples/example_executable`


> **Note:** This **FindCUGL.cmake** enables CMake to automate the usage of CUGL
library. For more information, check CMake tutorials or official documentation.
