set(ZEPHYR_PHYS_EOS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/eos.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ideal_gas.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mie_gruneisen.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/stiffened_gas.cpp
        PARENT_SCOPE)

install(FILES
        eos.h
        ideal_gas.h
        mie_gruneisen.h
        stiffened_gas.h
        DESTINATION
        include/${PROJECT_NAME}/phys/matter/eos)