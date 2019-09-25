cmake_minimum_required(VERSION 3.5)
project(oicugl VERSION 0.1.0)

#add_executable(bouncingball bouncingball.cpp)
add_executable(flyingplane flyingplane.cpp)
add_executable(gimbal gimbal.cpp)
add_executable(integration integration.cpp)

include_directories(${CMAKE_SOURCE_DIR}/include)

# Graphics APIs
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(CUGL REQUIRED)

if(GLUT_FOUND)
	#target_link_libraries(bouncingball ${GLUT_LIBRARIES})
	target_link_libraries(flyingplane ${GLUT_LIBRARIES})
	target_link_libraries(gimbal ${GLUT_LIBRARIES})
	target_link_libraries(integration ${GLUT_LIBRARIES})
endif()

if(OPENGL_FOUND)
	#target_link_libraries(bouncingball ${OPENGL_LIBRARIES})
	target_link_libraries(flyingplane ${OPENGL_LIBRARIES})
	target_link_libraries(gimbal ${OPENGL_LIBRARIES})
	target_link_libraries(integration ${OPENGL_LIBRARIES})
endif()

if(CUGL_LIBRARY_FOUND)
	#target_link_libraries(bouncingball ${CUGL_LIBRARY})
	target_link_libraries(flyingplane ${CUGL_LIBRARY})
	target_link_libraries(gimbal ${CUGL_LIBRARY})
	target_link_libraries(integration ${CUGL_LIBRARY})
endif()

# EOF
