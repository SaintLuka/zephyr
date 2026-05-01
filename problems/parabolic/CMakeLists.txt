list(APPEND PROBLEMS poisson)

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
    message(STATUS "Build parabolic problems:\n   ${PROBLEMS}")
endif()

foreach(problem IN LISTS PROBLEMS)
    add_executable(${problem} ${problem}.cpp)
    target_include_directories(${problem} PUBLIC zephyr)
    target_link_libraries(${problem} PUBLIC zephyr)
    install(TARGETS ${problem} DESTINATION problems/parabolic)
endforeach()