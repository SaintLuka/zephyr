list(APPEND TESTS dense)
#list(APPEND TESTS gradient)
#list(APPEND TESTS interface)
#list(APPEND TESTS random)
#list(APPEND TESTS templates)

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
    list(REMOVE_ITEM TESTS templates)
endif()

if (TESTS)
    message(STATUS "Build AMR tests:\n   ${TESTS}")
endif()

foreach(test IN LISTS TESTS)
    add_executable(${test} ${test}.cpp)
    target_include_directories(${test} PUBLIC zephyr)
    target_link_libraries(${test} PUBLIC zephyr)
    install(TARGETS ${test} DESTINATION tests/amr)
endforeach()
