set(ZEPHYR_MATH_FLUX_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/cir.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rusanov.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/godunov.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hllc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hll.cpp
        PARENT_SCOPE)

install(FILES
        num_flux.h
        cir.h
        rusanov.h
        godunov.h
        hllc.h
        hll.h
        DESTINATION
        include/${PROJECT_NAME}/math/cfd/fluxes)