#list(APPEND TESTS crp_interface)
#list(APPEND TESTS cuboid_clip)
#list(APPEND TESTS curves)
#list(APPEND TESTS mapping)
#list(APPEND TESTS poly_integrate)
#list(APPEND TESTS polygon_clip)
#list(APPEND TESTS polyhedra)
#list(APPEND TESTS splines)

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

# Can't compile targets without python
if (NOT ZEPHYR_PYTHON)
    list(REMOVE_ITEM TESTS crp_interface)
    list(REMOVE_ITEM TESTS curves)
    list(REMOVE_ITEM TESTS polygon_clip)
    list(REMOVE_ITEM TESTS splines)
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
