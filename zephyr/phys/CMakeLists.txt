add_subdirectory(eos)
add_subdirectory(tests)

set(ZEPHYR_PHYS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/fractions.cpp
        ${ZEPHYR_PHYS_EOS_SRC}
        ${ZEPHYR_PHYS_TESTS_SRC}
        PARENT_SCOPE)

install(FILES
        fractions.h
        literals.h
        DESTINATION
        include/${PROJECT_NAME}/physics)
