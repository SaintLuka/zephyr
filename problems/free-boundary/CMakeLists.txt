list(APPEND PROBLEMS fb_test_1D)

# Assuming all files in the directory are tests
if (ZEPHYR_PROBLEMS_ALL)
    file(GLOB PROBLEM_SOURCES *.cpp)

    set(PROBLEMS "")
    foreach(filename IN LISTS PROBLEM_SOURCES)
        # short name without path and extension
        get_filename_component(problem ${filename} NAME_WE)
        list(APPEND PROBLEMS ${problem})
    endforeach()
endif()

if (PROBLEMS)
    message(STATUS "Build free-boundary problems:\n   ${PROBLEMS}")
endif()

foreach(problem IN LISTS PROBLEMS)
    add_executable(${problem} ${problem}.cpp)
    target_include_directories(${problem} PUBLIC zephyr)
    target_link_libraries(${problem} PUBLIC zephyr)
    install(TARGETS ${problem} DESTINATION problems/free-boundary)
endforeach()
