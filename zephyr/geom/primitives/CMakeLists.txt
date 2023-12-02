set(ZEPHYR_GEOM_PRIMITIVES_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/bfaces.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/bvertices.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/amr_cell.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/amr_cell_check.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mov_cell.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/mov_cell_check.cxx
        PARENT_SCOPE
)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/base.h
        ${CMAKE_CURRENT_SOURCE_DIR}/adjacent.h
        ${CMAKE_CURRENT_SOURCE_DIR}/boundary.h
        ${CMAKE_CURRENT_SOURCE_DIR}/bface.h
        ${CMAKE_CURRENT_SOURCE_DIR}/bfaces.h
        ${CMAKE_CURRENT_SOURCE_DIR}/bvertices.h
        ${CMAKE_CURRENT_SOURCE_DIR}/amr_cell.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mov_node.h
        ${CMAKE_CURRENT_SOURCE_DIR}/bnodes.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mov_cell.h
        DESTINATION include/${PROJECT_NAME}/geom/primitives
)
