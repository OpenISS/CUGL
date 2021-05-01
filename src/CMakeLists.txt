#
AUX_SOURCE_DIRECTORY(. SOURCEFILES)
add_library(oicugl SHARED ${SOURCEFILES})

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
    		PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/include/cugl.h"
		)

# Graphics APIs
find_package(glew REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

target_link_libraries(oicugl PRIVATE GLEW::GLEW glm)

#
install(
	TARGETS oicugl
	ARCHIVE DESTINATION lib
    	LIBRARY DESTINATION lib
    	PUBLIC_HEADER DESTINATION include/cugl
    	)

# EOF
