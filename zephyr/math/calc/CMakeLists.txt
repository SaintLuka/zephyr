set(ZEPHYR_MATH_CALC_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/integrals.cpp
    PARENT_SCOPE)

install(FILES
        integrals.h
        DESTINATION
        include/${PROJECT_NAME}/math/calc)