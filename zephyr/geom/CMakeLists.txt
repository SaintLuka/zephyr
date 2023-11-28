add_subdirectory(generator)
add_subdirectory(primitives)

set(ZEPHYR_GEOM_SRC
        ${ZEPHYR_GEOM_GENERATOR_SRC}
        ${ZEPHYR_GEOM_PRIMITIVES_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/box.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/maps.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/polygon.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/grid.cpp
        PARENT_SCOPE
)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/vector.h
        ${CMAKE_CURRENT_SOURCE_DIR}/box.h
        ${CMAKE_CURRENT_SOURCE_DIR}/geom.h
        ${CMAKE_CURRENT_SOURCE_DIR}/maps.h
        ${CMAKE_CURRENT_SOURCE_DIR}/polygon.h
        ${CMAKE_CURRENT_SOURCE_DIR}/grid.h
        DESTINATION include/${PROJECT_NAME}/geom
)
