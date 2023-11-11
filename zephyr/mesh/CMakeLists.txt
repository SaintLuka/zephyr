add_subdirectory(euler)

set(ZEPHYR_MESH_SRC
    ${ZEPHYR_MESH_EULER_SRC}
    PARENT_SCOPE)

install(FILES
    storage.h
    threads.h
    DESTINATION
    include/${PROJECT_NAME}/mesh)
