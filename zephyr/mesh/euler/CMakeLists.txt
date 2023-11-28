set(ZEPHYR_MESH_EULER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/distributor.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/eu_face.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/eu_cell.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/eu_mesh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/eu_mesh_check.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/eu_mesh_refine.cpp
        PARENT_SCOPE)

install(FILES
    distributor.h
    eu_face.h
    eu_cell.h
    eu_mesh.h
    DESTINATION
    include/${PROJECT_NAME}/mesh/euler)
