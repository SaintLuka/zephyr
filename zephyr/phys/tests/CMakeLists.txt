set(ZEPHYR_PHYS_TESTS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/toro.cxx
        PARENT_SCOPE)

install(FILES
        shu_osher.h
        sod.h
        toro.h
        classic_test.h
        sedov.h
        mach.h
        RiemannTest2D.h
        quirck.h
        supersonic_flow_around_cylinder.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)