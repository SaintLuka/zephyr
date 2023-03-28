set(ZEPHYR_REFINER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/refiner.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/impl/check.cxx
    PARENT_SCOPE
)

install(FILES
    distributor.h
    refiner.h
    DESTINATION
    include/${PROJECT_NAME}/mesh/refiner
)