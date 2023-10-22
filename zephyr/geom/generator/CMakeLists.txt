set(ZEPHYR_GENERATOR_COLLECTION
	${CMAKE_CURRENT_SOURCE_DIR}/collection/plane_with_hole.cpp)

set(ZEPHYR_GEOM_GENERATOR_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/generator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/strip.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/rectangle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sector.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/cuboid.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/curve/curve.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/curve/plane.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/curve/circle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/curve/cubic.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bs_vertex.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/block.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/block_structured.cpp

	${ZEPHYR_GENERATOR_COLLECTION}
	PARENT_SCOPE
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/curve/circle.h
	${CMAKE_CURRENT_SOURCE_DIR}/curve/cubic.h
	${CMAKE_CURRENT_SOURCE_DIR}/curve/curve.h
        ${CMAKE_CURRENT_SOURCE_DIR}/curve/plane.h
	DESTINATION
	include/${PROJECT_NAME}/geom/generator/curve
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/generator.h
	${CMAKE_CURRENT_SOURCE_DIR}/strip.h
	${CMAKE_CURRENT_SOURCE_DIR}/rectangle.h
	${CMAKE_CURRENT_SOURCE_DIR}/sector.h
	${CMAKE_CURRENT_SOURCE_DIR}/cuboid.h
	${CMAKE_CURRENT_SOURCE_DIR}/bs_vertex.h
	${CMAKE_CURRENT_SOURCE_DIR}/block.h
	${CMAKE_CURRENT_SOURCE_DIR}/block_structured.h
	DESTINATION
	include/${PROJECT_NAME}/geom/generator
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/collection/plane_with_hole.h
	DESTINATION
	include/${PROJECT_NAME}/geom/generator/collection
	)
