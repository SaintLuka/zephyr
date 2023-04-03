set(ZEPHYR_MESH_REFINER_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/refiner.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/distributor.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/impl/check.cxx
    PARENT_SCOPE
)

install(FILES
    distributor.h
    refiner.h
    DESTINATION
    include/${PROJECT_NAME}/mesh/refiner
)