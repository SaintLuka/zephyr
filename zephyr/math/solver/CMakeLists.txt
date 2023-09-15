set(ZEPHYR_MATH_SOLVER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/riemann.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/transfer.cpp
        PARENT_SCOPE)

install(FILES
        convection.h
        riemann.h
        transfer.h
        DESTINATION
        include/${PROJECT_NAME}/math/solver)