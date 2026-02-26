#list(APPEND TESTS consistency)
#list(APPEND TESTS matter)
#list(APPEND TESTS mixture)
#list(APPEND TESTS mixture2)
#list(APPEND TESTS mixture3)

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
    list(REMOVE_ITEM TESTS mixture2)
    list(REMOVE_ITEM TESTS mixture3)
endif()

if (TESTS)
    message(STATUS "Build matter tests:\n   ${TESTS}")
endif()

foreach(test IN LISTS TESTS)
    add_executable(${test} ${test}.cpp)
    target_include_directories(${test} PUBLIC zephyr)
    target_link_libraries(${test} PUBLIC zephyr)
    install(TARGETS ${test} DESTINATION tests/eos)
endforeach()
