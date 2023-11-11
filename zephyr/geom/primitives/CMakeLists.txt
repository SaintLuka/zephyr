set(ZEPHYR_GEOM_PRIMITIVES_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/amr_faces.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/amr_cell.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/amr_cell_check.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_cell.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_cell_check.cxx
    PARENT_SCOPE
)

install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/base.h
    ${CMAKE_CURRENT_SOURCE_DIR}/adjacent.h
    ${CMAKE_CURRENT_SOURCE_DIR}/boundary.h
    ${CMAKE_CURRENT_SOURCE_DIR}/basic_face.h
    ${CMAKE_CURRENT_SOURCE_DIR}/amr_vertices.h
    ${CMAKE_CURRENT_SOURCE_DIR}/eu_face.h
    ${CMAKE_CURRENT_SOURCE_DIR}/amr_faces.h
    ${CMAKE_CURRENT_SOURCE_DIR}/amr_cell.h
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_node.h
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_nodes.h
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_face.h
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_faces.h
    ${CMAKE_CURRENT_SOURCE_DIR}/poly_cell.h
    DESTINATION include/${PROJECT_NAME}/geom/primitives
)
