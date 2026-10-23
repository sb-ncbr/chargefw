set(CHARGEFW_PYTHON_TESTS)
foreach(name IN ITEMS import molecule calculation gemmi output rdkit recipes)
    add_test(
        NAME test_chargefw_python_${name}
        COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" -m pytest -p no:cacheprovider -q
                "${CMAKE_CURRENT_SOURCE_DIR}/test_${name}.py"
    )
    list(APPEND CHARGEFW_PYTHON_TESTS test_chargefw_python_${name})
endforeach()

add_test(
    NAME test_chargefw_python_component_install
    COMMAND ${CMAKE_COMMAND}
            -DCHARGEFW_BINARY_DIR=${PROJECT_BINARY_DIR}
            -DCHARGEFW_EXPECTED_VERSION=${CHARGEFW_VERSION}
            -DCHARGEFW_PYTHON_EXECUTABLE=${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}
            -DCHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY=${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}
            -DCHARGEFW_TEST_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_python_component_install.cmake
)

set_tests_properties(${CHARGEFW_PYTHON_TESTS} PROPERTIES
    ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python;CHARGEFW_EXPECTED_VERSION=${CHARGEFW_VERSION}"
)
list(APPEND CHARGEFW_PYTHON_TESTS test_chargefw_python_component_install)
set_tests_properties(${CHARGEFW_PYTHON_TESTS} PROPERTIES FAIL_REGULAR_EXPRESSION "nanobind: leaked")

add_test(
    NAME test_chargefw_python_mypy
    COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" -m mypy --strict python/chargefw tests/python
)
set_tests_properties(test_chargefw_python_mypy PROPERTIES
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    ENVIRONMENT "MYPYPATH=${PROJECT_SOURCE_DIR}/python"
)
list(APPEND CHARGEFW_PYTHON_TESTS test_chargefw_python_mypy)

# A shared library loaded by an unsanitized Python interpreter cannot resolve Clang's sanitizer
# symbols unless the runtime is loaded into the process first. Native test executables link the runtime
# directly and do not need this adjustment.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
   (CHARGEFW_ENABLE_ASAN OR CHARGEFW_ENABLE_UBSAN))
    execute_process(
        COMMAND "${CMAKE_CXX_COMPILER}" --print-resource-dir
        OUTPUT_VARIABLE CHARGEFW_CLANG_RESOURCE_DIRECTORY
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )
    if(CHARGEFW_ENABLE_ASAN)
        set(CHARGEFW_SANITIZER_RUNTIME_NAME "libclang_rt.asan.so")
        set(CHARGEFW_SANITIZER_OPTIONS "ASAN_OPTIONS=detect_leaks=0")
    else()
        set(CHARGEFW_SANITIZER_RUNTIME_NAME "libclang_rt.ubsan_standalone.so")
        set(CHARGEFW_SANITIZER_OPTIONS "UBSAN_OPTIONS=print_stacktrace=1")
    endif()
    file(GLOB CHARGEFW_SANITIZER_RUNTIMES
         "${CHARGEFW_CLANG_RESOURCE_DIRECTORY}/lib/*/${CHARGEFW_SANITIZER_RUNTIME_NAME}")
    list(LENGTH CHARGEFW_SANITIZER_RUNTIMES CHARGEFW_SANITIZER_RUNTIME_COUNT)
    if(NOT CHARGEFW_SANITIZER_RUNTIME_COUNT EQUAL 1)
        message(FATAL_ERROR
            "Expected one ${CHARGEFW_SANITIZER_RUNTIME_NAME} below "
            "${CHARGEFW_CLANG_RESOURCE_DIRECTORY}, found ${CHARGEFW_SANITIZER_RUNTIME_COUNT}")
    endif()
    list(GET CHARGEFW_SANITIZER_RUNTIMES 0 CHARGEFW_SANITIZER_RUNTIME)
    set_property(TEST ${CHARGEFW_PYTHON_TESTS}
        APPEND PROPERTY ENVIRONMENT
        "LD_PRELOAD=${CHARGEFW_SANITIZER_RUNTIME}"
        "${CHARGEFW_SANITIZER_OPTIONS}"
    )
endif()
