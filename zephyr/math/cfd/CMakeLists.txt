add_subdirectory(flux)

set(ZEPHYR_MATH_CFD_SRC
        ${ZEPHYR_MATH_FLUX_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/models.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rotate.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/limiter.cxx
    PARENT_SCOPE)

install(FILES
        models.h
        rotate.h
        limiter.h
        face_extra.h
        gradient.h
        DESTINATION
        include/${PROJECT_NAME}/math/cfd)