set(ZEPHYR_IO_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/vtk_type.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/variable.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/variables.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/csv_file.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/vtu_file.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pvd_file.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/write_freq.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/vtk_type.h
        ${CMAKE_CURRENT_SOURCE_DIR}/variable.h
        ${CMAKE_CURRENT_SOURCE_DIR}/variables.h
        ${CMAKE_CURRENT_SOURCE_DIR}/csv_file.h
        ${CMAKE_CURRENT_SOURCE_DIR}/vtu_file.h
        ${CMAKE_CURRENT_SOURCE_DIR}/pvd_file.h
        ${CMAKE_CURRENT_SOURCE_DIR}/write_freq.h
        DESTINATION include/${PROJECT_NAME}/io)
