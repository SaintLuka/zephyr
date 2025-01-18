set(ZEPHYR_MATH_CALC_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/derivatives.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/integrals.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/weno.cpp
    PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/derivatives.h
        ${CMAKE_CURRENT_SOURCE_DIR}/integrals.h
        ${CMAKE_CURRENT_SOURCE_DIR}/weno.h
        DESTINATION
        include/${PROJECT_NAME}/math/calc)