# FindCUGL.cmake
# s_jashan@cs.concordia.ca

#
find_path(
	CUGL_INCLUDE_DIR
	NAMES "cugl/cugl.h"
	PATHS "$ENV{CUGL_INCLUDE}"
	)

find_library(
	CUGL_LIBRARY
	NAMES "oicugl"
	PATHS "$ENV{CUGL_REDIST}"
	)

set(CUGL_INCLUDE_DIR ${CUGL_INCLUDE_DIR})
set(CUGL_LIBRARY ${CUGL_LIBRARY})

if(CUGL_INCLUDE_DIR)
	message("CUGL_INCLUDE_DIR found at " ${CUGL_INCLUDE_DIR} )
	set (CUGL_INCLUDE_DIR_FOUND TRUE)
	include_directories(${CUGL_INCLUDE_DIR})
endif(CUGL_INCLUDE_DIR)

if(CUGL_LIBRARY)
	message("CUGL_LIBRARY found at " ${CUGL_LIBRARY} )
	set (CUGL_LIBRARY_FOUND TRUE)
	link_directories(${CUGL_LIBRARY})
endif(CUGL_LIBRARY)
#
