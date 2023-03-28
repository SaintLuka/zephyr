add_subdirectory(cfd)
#add_subdirectory(sle)

set(ZEPHYR_MATH_SRC
    ${ZEPHYR_MATH_CFD_SRC}
    PARENT_SCOPE)

install(FILES
    vectorization.h
    DESTINATION
    include/${PROJECT_NAME}/math)
