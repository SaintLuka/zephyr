list(APPEND PROBLEMS swe_drop)
list(APPEND PROBLEMS swe_tests)

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
    message(STATUS "Build shallow-water problems:\n   ${PROBLEMS}")
endif()

foreach(problem IN LISTS PROBLEMS)
    add_executable(${problem} ${problem}.cpp)
    target_include_directories(${problem} PUBLIC zephyr)
    target_link_libraries(${problem} PUBLIC zephyr)
    install(TARGETS ${problem} DESTINATION problems/single-fluid)
endforeach()
