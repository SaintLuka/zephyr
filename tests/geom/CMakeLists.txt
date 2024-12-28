#list(APPEND TESTS crp_interface)
#list(APPEND TESTS curves)
#list(APPEND TESTS game_of_life)
#list(APPEND TESTS maps_integral)
#list(APPEND TESTS mesh_generator)
#list(APPEND TESTS movable_mesh)
#list(APPEND TESTS polygon_clip)
#list(APPEND TESTS splines)
#list(APPEND TESTS volume_fraction)

# Assume all files in direcory are tests
if (ZEPHYR_TESTS_ALL)
    file(GLOB TEST_SOURCES *.cpp)

    set(TESTS "")
    foreach(filename IN LISTS TEST_SOURCES)
        # short name without path and extentsion
        get_filename_component(test ${filename} NAME_WE)
        list(APPEND TESTS ${test})
    endforeach()
endif()

if (TESTS)
    message(STATUS "Build geom tests:\n   ${TESTS}")
endif()

foreach(test IN LISTS TESTS)
    add_executable(${test} ${test}.cpp)
    target_include_directories(${test} PUBLIC zephyr)
    target_link_libraries(${test} PUBLIC zephyr)
    install(TARGETS ${test} DESTINATION tests/geom)
endforeach()
