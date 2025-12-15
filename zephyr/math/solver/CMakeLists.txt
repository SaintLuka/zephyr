set(ZEPHYR_MATH_SOLVER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sm_fluid.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mm_fluid.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/riemann.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/transfer.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.h
        ${CMAKE_CURRENT_SOURCE_DIR}/sm_fluid.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mm_fluid.h
        ${CMAKE_CURRENT_SOURCE_DIR}/riemann.h
        ${CMAKE_CURRENT_SOURCE_DIR}/transfer.h
        DESTINATION
        include/${PROJECT_NAME}/math/solver)