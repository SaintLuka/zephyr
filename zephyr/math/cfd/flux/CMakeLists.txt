set(ZEPHYR_MATH_FLUX_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/cir.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rusanov.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/godunov.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hllc.cpp
        PARENT_SCOPE)

install(FILES
        num_flux.h
        cir.h
        rusanov.h
        godunov.h
        hllc.h
        DESTINATION
        include/${PROJECT_NAME}/math/cfd/fluxes)