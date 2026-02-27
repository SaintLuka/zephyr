cmake_minimum_required(VERSION 3.14)
project(zephyr
    VERSION   0.1
    LANGUAGES CXX)

set(CMAKE_LIB_NAME zephyr)

# Флаги компиляции
set(CMAKE_CXX_STANDARD 20)
if (MSVC)
    set(BASIC_FLAGS "/W1")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${BASIC_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} ${BASIC_FLAGS} /Od")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${BASIC_FLAGS} /O2")
else()
    set(BASIC_FLAGS "-Wno-attributes -Wno-unused-result")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${BASIC_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} ${BASIC_FLAGS} -O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${BASIC_FLAGS} -O3")
endif()

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type: Release / Debug")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)

option(ZEPHYR_ASSERTS      "Enable asserts"     OFF)
option(ZEPHYR_DOXYGEN      "Create docs"        OFF)
option(ZEPHYR_EIGEN        "Use Eigen library"  ON)
option(ZEPHYR_MPI          "Use MPI functions"  OFF)
option(ZEPHYR_PROBLEMS     "Build problems"     ON)
option(ZEPHYR_PROBLEMS_ALL "Build ALL problems" OFF)
option(ZEPHYR_PYTHON       "Use python libs"    OFF)
option(ZEPHYR_TESTS        "Build tests"        ON)
option(ZEPHYR_TESTS_ALL    "Build ALL tests"    OFF)

set(THREADS_TYPE "STD" CACHE STRING "Choose thread style: STD, TBB")
set_property(CACHE THREADS_TYPE PROPERTY STRINGS STD TBB)
if (${THREADS_TYPE} STREQUAL TBB)
    message(STATUS "Use TBB")
    set(ZEPHYR_TBB ON)
else ()
    message(STATUS "Use standard threads")
    set(ZEPHYR_TBB OFF)
endif ()
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set(ZEPHYR_ASSERTS ON)
endif ()

if (ZEPHYR_TBB)
    find_package(TBB REQUIRED)
else ()
    find_package(Threads REQUIRED)
endif ()
find_package(Eigen3  REQUIRED)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost   COMPONENTS program_options REQUIRED)
if (ZEPHYR_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
    find_package(pybind11 REQUIRED)

    get_filename_component(PYTHON_EXE_DIR "${Python3_EXECUTABLE}" DIRECTORY)
    get_filename_component(PYTHONHOME "${PYTHON_EXE_DIR}" ABSOLUTE)

    if(WIN32)
        file(TO_CMAKE_PATH "${PYTHONHOME}" PYTHONHOME_CMAKE)
        string(REPLACE "/" "\\" PYTHONHOME_NATIVE "${PYTHONHOME_CMAKE}")
        string(REPLACE "\\" "\\\\" PYTHONHOME_CPP "${PYTHONHOME_NATIVE}")
        set(PYTHONHOME "${PYTHONHOME_CPP}")
    endif()

    message(STATUS "PYTHONHOME: ${PYTHONHOME}")
endif()

add_definitions(-D_USE_MATH_DEFINES)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/configuration.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/zephyr/configuration.h
)

add_subdirectory(zephyr)

# Тип компиляции библиотеки, потом сделаем изменяемым вместо хардкода
set(BUILD_SHARED_LIBS YES)

# Тип компиляции (SHARED/STATIC) пропущен, полагаемся на переменную BUILD_SHARED_LIBS
add_library(zephyr ${ZEPHYR_SOURCES})
if (ZEPHYR_TBB)
    target_link_libraries(zephyr PUBLIC TBB::tbb)
else ()
    target_link_libraries(zephyr PUBLIC ${CMAKE_THREAD_LIBS_INIT})
endif ()
target_link_libraries(zephyr PUBLIC ${Boost_LIBRARIES})
if (ZEPHYR_PYTHON)
    target_link_libraries(zephyr PUBLIC Python3::Python Python3::NumPy)
    target_link_libraries(zephyr PUBLIC pybind11::embed)
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

if(MSVC)
    set_target_properties(zephyr PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS TRUE
            CXX_VISIBILITY_PRESET hidden)

    #set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

    target_compile_definitions(zephyr PRIVATE ZEPHYR_EXPORTS)

    set_target_properties(zephyr PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )
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

if (${ZEPHYR_TESTS} OR ${ZEPHYR_TESTS_ALL})
    add_subdirectory(tests)
endif()

if (${ZEPHYR_PROBLEMS} OR ${ZEPHYR_PROBLEMS_ALL})
    add_subdirectory(problems)
endif()

# Включить опцию ZEPHYR_DOXYGEN, затем вызвать make docs
if (ZEPHYR_DOXYGEN)
    # Поиск Doxygen и Graphviz
    find_package(Doxygen REQUIRED)
    find_program(DOXYGEN_DOT_EXECUTABLE dot REQUIRED)

    if(NOT DOXYGEN_FOUND)
        message(WARNING "Doxygen not found! Documentation will not be generated.")
    else()
        # Настройка Doxygen
        set(DOXYGEN_INPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
        set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs")
        set(DOXYGEN_HTML_OUTPUT "${DOXYGEN_OUTPUT_DIR}/html")
        set(DOXYGEN_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")

        # Конфигурация Doxyfile
        configure_file(
                "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
                "${DOXYGEN_CONFIG_FILE}"
                @ONLY
        )

        # Добавление цели для генерации документации
        add_custom_target(docs
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_FILE}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating documentation with Doxygen"
                VERBATIM
        )
    endif()
endif()