set(ZEPHYR_PHYS_COND_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/conductivity.cpp
        PARENT_SCOPE)

install(FILES
        conductivity.h
        DESTINATION
        include/${PROJECT_NAME}/phys/matter/cond)