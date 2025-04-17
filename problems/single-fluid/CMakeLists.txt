#list(APPEND PROBLEMS clogged_shock_tube)
#list(APPEND PROBLEMS smf_flow_around)
#list(APPEND PROBLEMS smf_test_1D)
list(APPEND PROBLEMS smf_test_2D)
list(APPEND PROBLEMS smf_test_3D)
list(APPEND PROBLEMS smf_test_3D_soa)
#list(APPEND PROBLEMS smf_test_axial)

# Assume all files in direcory are tests
if (ZEPHYR_PROBLEMS_ALL)
    file(GLOB PROBLEM_SOURCES *.cpp)

    set(PROBLEMS "")
    foreach(filename IN LISTS PROBLEM_SOURCES)
        # short name without path and extentsion
        get_filename_component(problem ${filename} NAME_WE)
        list(APPEND PROBLEMS ${problem})
    endforeach()
endif()

if (PROBLEMS)
    message(STATUS "Build single-fluid problems:\n   ${PROBLEMS}")
endif()

foreach(problem IN LISTS PROBLEMS)
    add_executable(${problem} ${problem}.cpp)
    target_include_directories(${problem} PUBLIC zephyr)
    target_link_libraries(${problem} PUBLIC zephyr)
    install(TARGETS ${problem} DESTINATION problems/single-fluid)
endforeach()
