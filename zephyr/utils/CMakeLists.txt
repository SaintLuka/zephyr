set(ZEPHYR_UTILS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/threads.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/stopwatch.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/logger.cpp
    PARENT_SCOPE)

install(FILES
    error_list.h
    mpi.h
    threads.h
    stopwatch.h
    logger.h
    DESTINATION
    include/${PROJECT_NAME}/utils)
