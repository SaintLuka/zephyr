set(CLOGGED_SHOCK_TUBE 1)

set(SMF_TEST_1D        1)
set(SMF_TEST_2D        1)

#set(R_TESTER    1)
#set(RIEMANN     1)
#set(TESTS2D     1)


if (CLOGGED_SHOCK_TUBE)
    message(STATUS "Build clogged shock tube test")
    add_executable(clogged-shock-tube clogged_shock_tube.cpp)
    target_link_directories(clogged-shock-tube PUBLIC zephyr)
    target_link_libraries(clogged-shock-tube PUBLIC zephyr)
    install(TARGETS clogged-shock-tube DESTINATION problems/single-fluid)
endif()

if (SMF_TEST_1D)
    message(STATUS "Build 1D single fluid test")
    add_executable(smf-test-1D smf_test_1D.cpp)
    target_link_directories(smf-test-1D PUBLIC zephyr)
    target_link_libraries(smf-test-1D PUBLIC zephyr)
    install(TARGETS smf-test-1D DESTINATION problems/single-fluid)
endif()

if (SMF_TEST_2D)
    message(STATUS "Build 2D single fluid test")
    add_executable(smf-test-2D smf_test_2D.cpp)
    target_link_directories(smf-test-2D PUBLIC zephyr)
    target_link_libraries(smf-test-2D PUBLIC zephyr)
    install(TARGETS smf-test-2D DESTINATION problems/single-fluid)
endif()




if (RIEMANN)
    message(STATUS "Build Riemann test")
    add_executable(riemann riemann.cpp)
    target_include_directories(riemann PUBLIC zephyr)
    target_link_libraries(riemann PUBLIC zephyr)
    install(TARGETS riemann DESTINATION problems/single-fluid)
endif()

if (R_TESTER)
    message(STATUS "Build Riemann tester")
    add_executable(r_tester riemann_tester.cpp)
    target_include_directories(r_tester PUBLIC zephyr)
    target_link_libraries(r_tester PUBLIC zephyr)
    install(TARGETS r_tester DESTINATION problems/single-fluid)
endif()

if (TESTS2D)
    message(STATUS "2D tests")
    add_executable(tests2D smf_test_2D.cpp)
    target_include_directories(tests2D PUBLIC zephyr)
    target_link_libraries(tests2D PUBLIC zephyr)
    install(TARGETS tests2D DESTINATION problems/single-fluid)
endif()