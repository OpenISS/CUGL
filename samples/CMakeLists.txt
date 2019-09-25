cmake_minimum_required(VERSION 3.5)
project(oicugl VERSION 0.1.0)

add_executable(flyingplane flyingplane.cpp)

target_include_directories( flyingplane
		PUBLIC
		${CMAKE_SOURCE_DIR}/include
		)

# Graphics APIs
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(CUGL REQUIRED)

if(GLUT_FOUND)
	target_link_libraries(flyingplane ${GLUT_LIBRARIES})
endif()

if(OPENGL_FOUND)
	target_link_libraries(flyingplane ${OPENGL_LIBRARIES})
endif()

if(CUGL_LIBRARY_FOUND)
	target_include_directories(flyingplane PUBLIC ${CUGL_INCLUDE_DIR})
	target_link_libraries(flyingplane ${CUGL_LIBRARY})
endif()

# EOF
