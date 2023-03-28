set(ZEPHYR_PHYS_EOS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/eos.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ideal_gas.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/stiffened_gas.cpp
        PARENT_SCOPE)

install(FILES
        eos.h
        ideal_gas.h
        stiffened_gas.h
        types.h
        DESTINATION
        include/${PROJECT_NAME}/phys/eos)