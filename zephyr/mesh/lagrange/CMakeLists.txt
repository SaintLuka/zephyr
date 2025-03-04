set(ZEPHYR_MESH_LAGRANGE_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/la_face.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/la_cell.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/la_mesh.cpp
        PARENT_SCOPE)

install(FILES
        la_face.h
        la_cell.h
        la_mesh.h
        la_node.h
        mov_storage.h
        DESTINATION
        include/${PROJECT_NAME}/mesh/lagrange)
