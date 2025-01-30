set(ZEPHYR_PHYS_TESTS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/ivp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_1D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_2D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sedov_blast.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/ivp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/test_1D.h
        ${CMAKE_CURRENT_SOURCE_DIR}/test_2D.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests)