set(ZEPHYR_PHYS_TESTS_SRC
        PARENT_SCOPE)

install(FILES
        shu-osher.h
        sod.h
        toro.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)