set(ZEPHYR_IO_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/vtk_type.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/variable.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/variables.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/csv_file.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/vtu_file.cxx
        ${CMAKE_CURRENT_SOURCE_DIR}/pvd_file.cxx
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/vtk_type.h
        ${CMAKE_CURRENT_SOURCE_DIR}/filter.h
        ${CMAKE_CURRENT_SOURCE_DIR}/variable.h
        ${CMAKE_CURRENT_SOURCE_DIR}/variables.h
        ${CMAKE_CURRENT_SOURCE_DIR}/csv_file.h
        ${CMAKE_CURRENT_SOURCE_DIR}/vtu_file.h
        ${CMAKE_CURRENT_SOURCE_DIR}/pvd_file.h
        DESTINATION include/${PROJECT_NAME}/io)
