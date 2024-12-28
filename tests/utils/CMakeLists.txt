#list(APPEND TESTS balancing)
#list(APPEND TESTS decomposition)
#list(APPEND TESTS mpitest)
#list(APPEND TESTS range)
#list(APPEND TESTS stopwatch)
#list(APPEND TESTS threads)
#list(APPEND TESTS vectors)

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
    message(STATUS "Build utils tests:\n   ${TESTS}")
endif()

foreach(test IN LISTS TESTS)
    add_executable(${test} ${test}.cpp)
    target_include_directories(${test} PUBLIC zephyr)
    target_link_libraries(${test} PUBLIC zephyr)
    install(TARGETS ${test} DESTINATION tests/utils)
endforeach()
