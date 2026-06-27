set(ZEPHYR_GEOM_SURFACE_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/solid_body.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/solid_body_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/solid_body_3d.cpp
        PARENT_SCOPE
)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/solid_body.h
        ${CMAKE_CURRENT_SOURCE_DIR}/solid_body_2d.h
        ${CMAKE_CURRENT_SOURCE_DIR}/solid_body_3d.h
        DESTINATION include/${PROJECT_NAME}/geom/surface
)
