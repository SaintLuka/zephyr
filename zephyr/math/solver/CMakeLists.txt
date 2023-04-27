set(ZEPHYR_MATH_SOLVER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/riemann.cpp
        PARENT_SCOPE)

install(FILES
        convection.h
        riemann.h
        DESTINATION
        include/${PROJECT_NAME}/math/solver)