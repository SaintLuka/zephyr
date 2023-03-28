set(ZEPHYR_GENERATOR_COLLECTION
	${CMAKE_CURRENT_SOURCE_DIR}/collection/plane_with_hole.cxx)

set(ZEPHYR_MESH_GENERATOR_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/base.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/cuboid.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/sector.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/rectangle.cxx

	${CMAKE_CURRENT_SOURCE_DIR}/curve/curve.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/curve/line.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/curve/circle.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/curve/cubic.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/vertex.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/block.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/block_structured.cxx

	${ZEPHYR_GENERATOR_COLLECTION}
	PARENT_SCOPE
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/curve/circle.h
	${CMAKE_CURRENT_SOURCE_DIR}/curve/cubic.h
	${CMAKE_CURRENT_SOURCE_DIR}/curve/curve.h
	${CMAKE_CURRENT_SOURCE_DIR}/curve/line.h
	DESTINATION
	include/${PROJECT_NAME}/mesh/generator/curve
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/base.h
	${CMAKE_CURRENT_SOURCE_DIR}/block.h
	${CMAKE_CURRENT_SOURCE_DIR}/block_structured.h
	${CMAKE_CURRENT_SOURCE_DIR}/cuboid.h
	${CMAKE_CURRENT_SOURCE_DIR}/rectangle.h
	${CMAKE_CURRENT_SOURCE_DIR}/sector.h
	${CMAKE_CURRENT_SOURCE_DIR}/vertex.h
	DESTINATION
	include/${PROJECT_NAME}/mesh/generator
)

install(FILES
	collection/plane_with_hole.h
	DESTINATION
	include/${PROJECT_NAME}/mesh/generator/collection
	)
