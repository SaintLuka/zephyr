add_subdirectory(decomp)
add_subdirectory(euler)

set(ZEPHYR_MESH_SRC
        ${ZEPHYR_MESH_EULER_SRC}
        ${ZEPHYR_MESH_DECOMP_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/storage.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/index.h
        ${CMAKE_CURRENT_SOURCE_DIR}/side.h
        ${CMAKE_CURRENT_SOURCE_DIR}/storage.h
        DESTINATION
        include/${PROJECT_NAME}/mesh)
