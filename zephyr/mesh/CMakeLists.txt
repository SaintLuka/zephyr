add_subdirectory(euler)
add_subdirectory(decomp)
add_subdirectory(primitives)

set(ZEPHYR_MESH_SRC
    ${ZEPHYR_MESH_EULER_SRC}
    ${ZEPHYR_MESH_DECOMP_SOURCES}
    ${ZEPHYR_MESH_PRIMITIVES_SRC}
    PARENT_SCOPE)

install(FILES
        storage.h
        DESTINATION
        include/${PROJECT_NAME}/mesh)
