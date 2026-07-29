set(ZEPHYR_SWE_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/dam_break.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/step.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/thacker.cpp
        PARENT_SCOPE)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/dam_break.h
        ${CMAKE_CURRENT_SOURCE_DIR}/step.h
        ${CMAKE_CURRENT_SOURCE_DIR}/sw_test.h
        ${CMAKE_CURRENT_SOURCE_DIR}/thacker.h
        DESTINATION
        include/${PROJECT_NAME}/phys/tests/swe)