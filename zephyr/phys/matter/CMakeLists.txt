add_subdirectory(cond)
add_subdirectory(eos)
add_subdirectory(plast)
add_subdirectory(visc)

set(ZEPHYR_PHYS_MATTER_SRC
        ${ZEPHYR_PHYS_COND_SRC}
        ${ZEPHYR_PHYS_EOS_SRC}
        ${ZEPHYR_PHYS_PLAST_SRC}
        ${ZEPHYR_PHYS_VISC_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/material.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/materials.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mixture_pt.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/material.h
        ${CMAKE_CURRENT_SOURCE_DIR}/materials.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mixture_pt.h
        DESTINATION
        include/${PROJECT_NAME}/phys/matter)