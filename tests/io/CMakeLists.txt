list(APPEND TESTS vtu_euler)
list(APPEND TESTS vtu_nodes)
list(APPEND TESTS vtu_polyhedra)
list(APPEND TESTS vtu_unique)
list(APPEND TESTS write_funcs)

# Assume all files in directory are tests
if (ZEPHYR_TESTS_ALL)
    file(GLOB TEST_SOURCES *.cpp)

    set(TESTS "")
    foreach(filename IN LISTS TEST_SOURCES)
        # short name without path and extension
        get_filename_component(test ${filename} NAME_WE)
        list(APPEND TESTS ${test})
    endforeach()
endif()

if (TESTS)
    message(STATUS "Build io tests:\n   ${TESTS}")
endif()

foreach(test IN LISTS TESTS)
    add_executable(${test} ${test}.cpp)
    target_include_directories(${test} PUBLIC zephyr)
    target_link_libraries(${test} PUBLIC zephyr)
    install(TARGETS ${test} DESTINATION tests/io)
endforeach()
