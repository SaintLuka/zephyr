cmake_minimum_required(VERSION 3.14)
project(zephyr
    VERSION   0.1
    LANGUAGES CXX)

set(CMAKE_LIB_NAME zephyr)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG   "-Wno-attributes -Wno-unused-result -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Wno-attributes -Wno-unused-result -O3")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(ZEPHYR_DOXYGEN      "Create docs"        OFF)
option(ZEPHYR_EIGEN        "Use Eigen library"  ON)
option(ZEPHYR_MPI          "Use MPI functions"  OFF)
option(ZEPHYR_PROBLEMS     "Build problems"     ON)
option(ZEPHYR_PROBLEMS_ALL "Build ALL problems" OFF)
option(ZEPHYR_PYTHON       "Use python libs"    ON)
option(ZEPHYR_TESTS        "Build tests"        ON)
option(ZEPHYR_TESTS_ALL    "Build ALL tests"    OFF)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/configuration.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/zephyr/configuration.h
)

find_package(Threads REQUIRED)
find_package(Eigen3  REQUIRED)
find_package(Boost   COMPONENTS program_options REQUIRED)
if (ZEPHYR_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
endif()

add_subdirectory(zephyr)

# Тип компиляции библиотеки, потом сделаем изменяемым вместо хардкода
set(BUILD_SHARED_LIBS YES)

# Тип компиляции (SHARED/STATIC) пропущен, полагаемся на переменную BUILD_SHARED_LIBS
add_library(zephyr ${ZEPHYR_SOURCES})
target_link_libraries(zephyr PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(zephyr PUBLIC ${Boost_LIBRARIES})
if (ZEPHYR_PYTHON)
    target_link_libraries(zephyr PUBLIC Python3::Python Python3::NumPy)
endif()
target_link_libraries(zephyr PUBLIC Eigen3::Eigen)

target_include_directories(zephyr PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${Boost_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    $<INSTALL_INTERFACE:${Boost_INCLUDE_DIR}>)

if (ZEPHYR_MPI)
    #find_package(MPI REQUIRED COMPONENTS C)
    #target_link_libraries(zephyr PUBLIC MPI::MPI_C)

    find_package(MPI REQUIRED)
    target_include_directories(zephyr PUBLIC ${MPI_INCLUDE_PATH})
    target_link_libraries(zephyr PUBLIC ${MPI_LIBRARIES})
endif()
install(TARGETS zephyr DESTINATION lib EXPORT zephyr-targets)

install(EXPORT zephyr-targets DESTINATION share/cmake/zephyr)

include(CMakePackageConfigHelpers)
# gen_eu the config file that is includes the exports
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

if (ZEPHYR_TESTS OR ZEPHYR_TESTS_ALL)
    add_subdirectory(tests)
endif()

if (ZEPHYR_PROBLEMS OR ZEPHYR_PROBLEMS_ALL)
    add_subdirectory(problems)
endif()

add_subdirectory(python)

if (ZEPHYR_DOXYGEN)
    find_package(Doxygen REQUIRED dot)
    doxygen_add_docs(docs ${PROJECT_SOURCE_DIR}/zephyr)
endif()