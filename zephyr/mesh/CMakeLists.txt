add_subdirectory(euler)
add_subdirectory(lagrange)
add_subdirectory(decomp)

set(ZEPHYR_MESH_SRC
    ${ZEPHYR_MESH_EULER_SRC}
    ${ZEPHYR_DECOMPOSITION_SOURCES}
    ${ZEPHYR_MESH_LAGRANGE_SRC}
    PARENT_SCOPE)

install(FILES
        storage.h
        face.h
        cell.h
        mesh.h
        DESTINATION
        include/${PROJECT_NAME}/mesh)
