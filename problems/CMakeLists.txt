add_subdirectory(single-fluid)
#add_subdirectory(multi-fluid)

#list(APPEND PROBLEMS advection)
#list(APPEND PROBLEMS convection)
#list(APPEND PROBLEMS membrane)
#list(APPEND PROBLEMS shu_osher)
#list(APPEND PROBLEMS transfer)
#list(APPEND PROBLEMS wave)

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
    message(STATUS "Build some problems:\n   ${PROBLEMS}")
endif()

foreach(problem IN LISTS PROBLEMS)
    add_executable(${problem} ${problem}.cpp
            wave.cpp)
    target_include_directories(${problem} PUBLIC zephyr)
    target_link_libraries(${problem} PUBLIC zephyr)
    install(TARGETS ${problem} DESTINATION problems)
endforeach()
