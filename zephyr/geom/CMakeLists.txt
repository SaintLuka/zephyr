set(ZEPHYR_GEOM_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/cell.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/cell_check.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/box.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/geom.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/maps.cxx
        PARENT_SCOPE
)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/vector.h
        ${CMAKE_CURRENT_SOURCE_DIR}/base.h
        ${CMAKE_CURRENT_SOURCE_DIR}/side.h
        ${CMAKE_CURRENT_SOURCE_DIR}/box.h
        ${CMAKE_CURRENT_SOURCE_DIR}/cell.h
        ${CMAKE_CURRENT_SOURCE_DIR}/face.h
        ${CMAKE_CURRENT_SOURCE_DIR}/faces.h
        ${CMAKE_CURRENT_SOURCE_DIR}/vertices.h
        ${CMAKE_CURRENT_SOURCE_DIR}/maps.h
        ${CMAKE_CURRENT_SOURCE_DIR}/geom.h
        DESTINATION include/${PROJECT_NAME}/geom
)
