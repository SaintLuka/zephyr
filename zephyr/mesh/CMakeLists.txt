add_subdirectory(generator)

set(ZEPHYR_MESH_SRC
    ${ZEPHYR_MESH_GENERATOR_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/storage.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/face.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/cell.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/range.cxx
    PARENT_SCOPE)


install(FILES
    mesh.h
    storage.h
    face.h
    cell.h
    range.h
    threads.h
    DESTINATION
    include/${PROJECT_NAME}/mesh)
