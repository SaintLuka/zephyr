#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "zephyr" for configuration "Debug"
set_property(TARGET zephyr APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(zephyr PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libzephyr.so"
  IMPORTED_SONAME_DEBUG "libzephyr.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS zephyr )
list(APPEND _IMPORT_CHECK_FILES_FOR_zephyr "${_IMPORT_PREFIX}/lib/libzephyr.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
