set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

aux_source_directory(. SRC_LIST)

add_executable(pnmatrix_test ${SRC_LIST})

add_test(NAME matrix_test COMMAND pnmatrix_test)
