set(ZEPHYR_MATH_SLE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/tridiagonal.cpp
    PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/tridiagonal.h
        DESTINATION
        include/${PROJECT_NAME}/math/sle)