set(ZEPHYR_MATH_SOLVER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/convection.cxx
        PARENT_SCOPE)

install(FILES
        convection.h
        DESTINATION
        include/${PROJECT_NAME}/math/solver)