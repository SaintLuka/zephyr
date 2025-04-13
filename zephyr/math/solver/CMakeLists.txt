set(ZEPHYR_MATH_SOLVER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sm_fluid.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sm_fluid_soa.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mm_fluid.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/riemann.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/transfer.cpp
        PARENT_SCOPE)

install(FILES
        convection.h
        sm_fluid.h
        sm_fluid_soa.h
        mm_fluid.h
        riemann.h
        transfer.h
        DESTINATION
        include/${PROJECT_NAME}/math/solver)