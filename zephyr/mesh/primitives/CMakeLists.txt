set(ZEPHYR_MESH_PRIMITIVES_SRC
        PARENT_SCOPE
)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/decomposition.h
        ${CMAKE_CURRENT_SOURCE_DIR}/side.h
        DESTINATION include/${PROJECT_NAME}/mesh/primitives
)
