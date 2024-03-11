set(ZEPHYR_PYTHON_FILES
        csv_file.py
        plot1D.py)

file(COPY
        ${ZEPHYR_PYTHON_FILES}
        DESTINATION
        ${CMAKE_BINARY_DIR}/python)

install(FILES
        ${ZEPHYR_PYTHON_FILES}
        DESTINATION python)