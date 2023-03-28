set(ZEPHYR_MATH_CFD_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/cir.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/exact.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/models.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rotate.cpp
        PARENT_SCOPE)

install(FILES
        cir.h
        exact.h
        fluxes.h
        hll.h
        models.h
        rotate.h
        DESTINATION
        include/${PROJECT_NAME}/math/cfd)