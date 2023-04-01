project(zephyr C CXX)
cmake_minimum_required(VERSION 3.5)

set(CMAKE_LIB_NAME zephyr)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG   "-Wno-attributes -Wno-unused-result -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Wno-attributes -Wno-unused-result -O3")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(ZEPHYR_ENABLE_MPI      "Use MPI functions" ON)
option(ZEPHYR_ENABLE_YAML     "YAML support"      OFF)
option(ZEPHYR_ENABLE_TESTS    "Build tests"       ON)
option(ZEPHYR_ENABLE_PROBLEMS "Build problems"    OFF)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/configuration.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/zephyr/configuration.h
)

find_package(Threads REQUIRED)
find_package(Boost COMPONENTS system filesystem)

if (ZEPHYR_ENABLE_YAML)
    find_package(yaml-cpp REQUIRED)
endif()

add_subdirectory(zephyr)

message(STATUS "zephyr: build shared")

add_library(zephyr SHARED ${ZEPHYR_SOURCES})
target_link_libraries(zephyr PUBLIC ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(zephyr PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/libs
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${Boost_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    $<INSTALL_INTERFACE:${Boost_INCLUDE_DIR}>)

if (ZEPHYR_ENABLE_MPI)
    #find_package(MPI REQUIRED COMPONENTS C)
    #target_link_libraries(zephyr PUBLIC MPI::MPI_C)

    find_package(MPI REQUIRED)
    target_include_directories(zephyr PUBLIC ${MPI_INCLUDE_PATH})
    target_link_libraries(zephyr PUBLIC ${MPI_LIBRARIES})
endif()
if (ZEPHYR_ENABLE_YAML)
    target_link_libraries(zephyr PUBLIC yaml-cpp)
endif()
install(TARGETS zephyr DESTINATION lib EXPORT zephyr-targets)

install(EXPORT zephyr-targets DESTINATION share/cmake/zephyr)
install(DIRECTORY libs/Eigen DESTINATION include)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/zephyr-config.cmake"
    INSTALL_DESTINATION "share/cmake/zephyr"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

# install the configuration file
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/zephyr-config.cmake
    DESTINATION share/cmake/zephyr
    )

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/zephyr/configuration.h
    DESTINATION include/${PROJECT_NAME}
    )

if (ZEPHYR_ENABLE_TESTS)
    add_subdirectory(tests)
endif()

if (ZEPHYR_ENABLE_PROBLEMS)
    add_subdirectory(problems)
endif()
