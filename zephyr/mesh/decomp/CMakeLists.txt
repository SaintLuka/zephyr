set(ZEPHYR_DECOMPOSITION_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/decomposition.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/orb/transform.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/orb/blocks.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/orb/ORB.cpp
	PARENT_SCOPE
)

install(FILES
		decomposition.h
		ORB.h
		VD3.h
		morton.h
		hilbert.h
		DESTINATION include/${PROJECT_NAME}/mesh/decomp
		)

install(FILES
		orb/blocks.h
		orb/transform.h
	DESTINATION include/${PROJECT_NAME}/network/decomp/ORB
)