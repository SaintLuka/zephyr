set(ZEPHYR_PHYS_TESTS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/toro.cxx
        PARENT_SCOPE)

install(FILES
        shu_osher.h
        sod.h
        toro.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)