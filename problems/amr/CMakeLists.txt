set(INTERFACE 1)

if (INTERFACE)
    message(STATUS "Build AMR interface test")
    add_executable(interface interface.cxx)
    target_include_directories(interface PUBLIC zephyr)
    target_link_libraries(interface PUBLIC zephyr)
    install(TARGETS interface DESTINATION problems/amr)
endif()