set(ZEPHYR_PHYS_TESTS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/toro.cxx
        PARENT_SCOPE)

install(FILES
        shu_osher.h
        sod.h
        toro.h
        classic_test.h
        blast_wave.h
        shock-in-a-box.h
        mach.h
        RiemannTest2D.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)