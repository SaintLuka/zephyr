#list(APPEND TESTS axial)
#list(APPEND TESTS balancing)
#list(APPEND TESTS distributed-amr)
#list(APPEND TESTS distributed-conv)
#list(APPEND TESTS distributed-conv-amr)
#list(APPEND TESTS game_of_life)
#list(APPEND TESTS mesh_generator)
#list(APPEND TESTS storage)
#list(APPEND TESTS voronoi)

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

# Can't compile targets without python
if (NOT ZEPHYR_PYTHON)
    list(REMOVE_ITEM TESTS voronoi)
endif()

if (TESTS)
    message(STATUS "Build utils tests:\n   ${TESTS}")
endif()

foreach(test IN LISTS TESTS)
    add_executable(${test} ${test}.cpp)
    target_include_directories(${test} PUBLIC zephyr)
    target_link_libraries(${test} PUBLIC zephyr)
    install(TARGETS ${test} DESTINATION tests/utils)
endforeach()
