set(ZEPHYR_GEOM_PRIMITIVES_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/line.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quad.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cube.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/triangle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polygon.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/polyhedron.cpp
    PARENT_SCOPE
)

install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/line.h
    ${CMAKE_CURRENT_SOURCE_DIR}/quad.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cube.h
    ${CMAKE_CURRENT_SOURCE_DIR}/triangle.h
    ${CMAKE_CURRENT_SOURCE_DIR}/polygon.h
    ${CMAKE_CURRENT_SOURCE_DIR}/polyhedron.h
    DESTINATION include/${PROJECT_NAME}/geom/primitives
)
