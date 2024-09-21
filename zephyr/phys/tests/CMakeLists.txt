set(ZEPHYR_PHYS_TESTS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/toro.cpp
        PARENT_SCOPE)

install(FILES
        shu_osher.h
        sod.h
        toro.h
        test_1D.h
        sedov.h
    shock_wave.h
        test_2D.h
    shock_wave.h
        rarefied_water.h
        richtmyer_meshkov.h
        supersonic_flow_around_cylinder.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)