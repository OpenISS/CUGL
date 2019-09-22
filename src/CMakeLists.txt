cmake_minimum_required(VERSION 3.12)
project(oicugl VERSION 0.1.0)

#
AUX_SOURCE_DIRECTORY(. SOURCEFILES)
add_library(oicugl ${SOURCEFILES})

target_include_directories( oicugl
		PUBLIC
		${CMAKE_SOURCE_DIR}/include
		)

set_target_properties(
		oicugl
		PROPERTIES
		ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
		LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
		)

# Graphics APIs
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)

if(GLUT_FOUND)
	target_link_libraries(oicugl ${GLUT_LIBRARIES})
endif()

if(OPENGL_FOUND)
	target_link_libraries(oicugl ${OPENGL_LIBRARIES})
endif()

# EOF