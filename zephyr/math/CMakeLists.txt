add_subdirectory(calc)
add_subdirectory(cfd)
add_subdirectory(sle)
add_subdirectory(solver)

set(ZEPHYR_MATH_SRC
    ${ZEPHYR_MATH_CALC_SRC}
    ${ZEPHYR_MATH_CFD_SRC}
    ${ZEPHYR_MATH_SLE_SRC}
    ${ZEPHYR_MATH_SOLVER_SRC}
    PARENT_SCOPE)

install(FILES
    vectorization.h
    DESTINATION
    include/${PROJECT_NAME}/math)
