add_subdirectory(decomp)
add_subdirectory(raw)

set(ZEPHYR_MESH_SRC
        ${ZEPHYR_MESH_RAW_SRC}
        ${ZEPHYR_MESH_DECOMP_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/storage.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/node.h
        ${CMAKE_CURRENT_SOURCE_DIR}/cell.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mesh.h
        ${CMAKE_CURRENT_SOURCE_DIR}/index.h
        ${CMAKE_CURRENT_SOURCE_DIR}/memory.h
        ${CMAKE_CURRENT_SOURCE_DIR}/storage.h
        DESTINATION
        include/${PROJECT_NAME}/mesh)
