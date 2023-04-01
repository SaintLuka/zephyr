add_subdirectory(geom)
add_subdirectory(io)
add_subdirectory(mesh)
add_subdirectory(math)
add_subdirectory(phys)
add_subdirectory(utils)

set(ZEPHYR_SOURCES
    ${ZEPHYR_GEOM_SRC}
    ${ZEPHYR_IO_SRC}
    ${ZEPHYR_MESH_SRC}
    ${ZEPHYR_MATH_SRC}
    ${ZEPHYR_PHYS_SRC}
    ${ZEPHYR_UTILS_SRC}
    PARENT_SCOPE
    )
