add_subdirectory(curves)
add_subdirectory(generator)

set(ZEPHYR_GEOM_SRC
    ${ZEPHYR_GEOM_CURVES_SRC}
    ${ZEPHYR_GEOM_GENERATOR_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/box.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/line.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quad.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cube.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sections.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/triangle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polygon.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polyhedron.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/intersection.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interface_recovery.cpp
    PARENT_SCOPE
)

install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary.h
    ${CMAKE_CURRENT_SOURCE_DIR}/box.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cell_type.h
    ${CMAKE_CURRENT_SOURCE_DIR}/line.h
    ${CMAKE_CURRENT_SOURCE_DIR}/quad.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cube.h
    ${CMAKE_CURRENT_SOURCE_DIR}/sections.h
    ${CMAKE_CURRENT_SOURCE_DIR}/triangle.h
    ${CMAKE_CURRENT_SOURCE_DIR}/polygon.h
    ${CMAKE_CURRENT_SOURCE_DIR}/polyhedron.h
    ${CMAKE_CURRENT_SOURCE_DIR}/geom.h
    ${CMAKE_CURRENT_SOURCE_DIR}/grid.h
    ${CMAKE_CURRENT_SOURCE_DIR}/intersection.h
    ${CMAKE_CURRENT_SOURCE_DIR}/interface_recovery.h
    ${CMAKE_CURRENT_SOURCE_DIR}/vector.h
    DESTINATION include/${PROJECT_NAME}/geom
)
