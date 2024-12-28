add_subdirectory(matter)
add_subdirectory(tests)

set(ZEPHYR_PHYS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/fractions.cpp
        ${ZEPHYR_PHYS_MATTER_SRC}
        ${ZEPHYR_PHYS_TESTS_SRC}
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/fractions.h
        ${CMAKE_CURRENT_SOURCE_DIR}/literals.h
        DESTINATION
        include/${PROJECT_NAME}/phys)
