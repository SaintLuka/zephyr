set(ZEPHYR_MATH_SOLVER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/riemann.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SmFluid.cpp
        PARENT_SCOPE)

install(FILES
        convection.h
        riemann.h
        SmFluid.h
        DESTINATION
        include/${PROJECT_NAME}/math/solver)