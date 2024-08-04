set(ZEPHYR_PHYS_TESTS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/toro.cpp
        PARENT_SCOPE)

install(FILES
        shu_osher.h
        sod.h
        toro.h
        test_1D.h
        sedov.h
        quirck.h
        test_2D.h
        quirck.h
        rarefied_water.h
        supersonic_flow_around_cylinder.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)