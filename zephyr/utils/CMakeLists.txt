set(ZEPHYR_UTILS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/threads.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/stopwatch.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/logger.cpp
    PARENT_SCOPE)

install(FILES
    error_list.h
    matplotlib.h
    mpi.h
    threads.h
    range.h
    stopwatch.h
    logger.h
    DESTINATION
    include/${PROJECT_NAME}/utils)
