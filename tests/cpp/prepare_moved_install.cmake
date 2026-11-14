include("${CMAKE_CURRENT_LIST_DIR}/support/test_cli.cmake")

if(NOT DEFINED CHARGEFW_BUILD_DIR OR NOT DEFINED CHARGEFW_INSTALL_PREFIX)
    message(FATAL_ERROR "Moved-install test setup requires a build directory and prefix")
endif()

set(install_source_prefix "${CHARGEFW_INSTALL_PREFIX}_source")
file(REMOVE_RECURSE "${install_source_prefix}" "${CHARGEFW_INSTALL_PREFIX}")

run_checked("ChargeFW installation"
            "${CMAKE_COMMAND}" --install "${CHARGEFW_BUILD_DIR}" --prefix "${install_source_prefix}")

file(RENAME "${install_source_prefix}" "${CHARGEFW_INSTALL_PREFIX}")
