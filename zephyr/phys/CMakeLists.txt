add_subdirectory(eos)
add_subdirectory(tests)

set(ZEPHYR_PHYS_SRC
        ${ZEPHYR_PHYS_EOS_SRC}
        ${ZEPHYR_PHYS_TESTS_SRC}
        PARENT_SCOPE)

install(FILES
        literals.h
        DESTINATION
        include/${PROJECT_NAME}/physics)
