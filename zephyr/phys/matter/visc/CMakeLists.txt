set(ZEPHYR_PHYS_VISC_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/viscosity.cpp
        PARENT_SCOPE)

install(FILES
        viscosity.h
        DESTINATION
        include/${PROJECT_NAME}/phys/matter/visc)