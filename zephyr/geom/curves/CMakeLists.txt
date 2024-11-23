set(ZEPHYR_GEOM_CURVES_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/bezier.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cubic_spline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/interpolant.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lagrange.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/linear_spline.cpp
        PARENT_SCOPE
)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/bezier.h
        ${CMAKE_CURRENT_SOURCE_DIR}/cubic_spline.h
        ${CMAKE_CURRENT_SOURCE_DIR}/interpolant.h
        ${CMAKE_CURRENT_SOURCE_DIR}/lagrange.h
        ${CMAKE_CURRENT_SOURCE_DIR}/linear_spline.h
        DESTINATION include/${PROJECT_NAME}/geom/curves
)
