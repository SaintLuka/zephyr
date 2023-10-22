set(ZEPHYR_PHYS_EOS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/eos.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ideal_gas.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/stiffened_gas.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mie_gruneisen.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/materials.cpp
        PARENT_SCOPE)

install(FILES
    types.h
    eos.h
    ideal_gas.h
    stiffened_gas.h
    mie_gruneisen.h
    materials.h
    DESTINATION
    include/${PROJECT_NAME}/phys/eos)