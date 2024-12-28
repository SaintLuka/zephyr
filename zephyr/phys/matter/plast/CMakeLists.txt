set(ZEPHYR_PHYS_PLAST_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/plasticity.cpp
        PARENT_SCOPE)

install(FILES
        plasticity.h
        DESTINATION
        include/${PROJECT_NAME}/phys/matter/plast)