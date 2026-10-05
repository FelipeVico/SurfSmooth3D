file(MAKE_DIRECTORY "${OUTPUT_DIR}")
execute_process(COMMAND "${RETENTION}" "${OPERATION}" "${REQUEST}"
  "${HELPER_REQUEST}" "${OUTPUT_DIR}/retained"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
file(WRITE "${OUTPUT_DIR}/result.json" "${output}")
file(WRITE "${OUTPUT_DIR}/stderr.txt" "${error}")
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Retained-context operation failed: ${error}")
endif()
string(JSON operations_equal GET "${output}" operations_byte_identical)
string(JSON measurements_equal GET "${output}" measurements_equal)
string(JSON helpers GET "${output}" helpers)
if(NOT operations_equal OR NOT measurements_equal OR helpers LESS 1)
  message(FATAL_ERROR "Retained CAD results changed, or fragmentation was not exercised: ${output}")
endif()
# Gmsh can add p-curve representations to a retained BREP during fragmentation.
# The existing test records those bytes; exact operation/measurement results
# are the established invariant, rather than byte-identical BREP serialization.
