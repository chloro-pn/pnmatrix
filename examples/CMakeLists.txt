set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

aux_source_directory(. SRC_LIST)

add_executable(pnmatrix_example ${SRC_LIST})
