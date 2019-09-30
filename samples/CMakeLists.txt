cmake_minimum_required(VERSION 3.12)
project(oicugl VERSION 0.1.0)

add_executable(sample sample.cpp)

target_include_directories( sample
		PUBLIC
		${CMAKE_SOURCE_DIR}/include
		)

# Graphics APIs
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(CUGL REQUIRED)

if(GLUT_FOUND)
	target_link_libraries(sample ${GLUT_LIBRARIES})
endif()

if(OPENGL_FOUND)
	target_link_libraries(sample ${OPENGL_LIBRARIES})
endif()

if(CUGL_LIBRARY_FOUND)
	target_include_directories(sample PUBLIC ${CUGL_INCLUDE_DIR})
	target_link_libraries(sample ${CUGL_LIBRARY})
endif()

# EOF