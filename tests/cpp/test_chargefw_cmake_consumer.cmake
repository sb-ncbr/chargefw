include("${CMAKE_CURRENT_LIST_DIR}/support/test_cli.cmake")

# The installation comes from the test_chargefw_install fixture and is only read here.
if(NOT DEFINED CHARGEFW_INSTALL_PREFIX OR NOT DEFINED CHARGEFW_PARENT_CACHE)
    message(FATAL_ERROR
            "Downstream consumer test requires the installation prefix and parent settings")
endif()

expect_files_exist("Installed third-party notices"
                   "${CHARGEFW_INSTALL_PREFIX}/${CHARGEFW_INSTALL_DOCDIR}/licenses/README.md")

set(consumer_build_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_build")
file(REMOVE_RECURSE "${consumer_build_directory}")

run_checked("Downstream ChargeFW consumer configuration"
            "${CMAKE_COMMAND}"
            -S "${CMAKE_CURRENT_LIST_DIR}/consumer"
            -B "${consumer_build_directory}"
            -C "${CHARGEFW_PARENT_CACHE}"
            -DCMAKE_PREFIX_PATH=${CHARGEFW_INSTALL_PREFIX}
            -DCMAKE_DISABLE_FIND_PACKAGE_nlohmann_json=ON
            -DCMAKE_DISABLE_FIND_PACKAGE_gemmi=ON)
run_checked("Downstream ChargeFW consumer build"
            "${CMAKE_COMMAND}" --build "${consumer_build_directory}")
run_checked("Downstream ChargeFW consumer" "${consumer_build_directory}/chargefw_consumer"
            "${CMAKE_CURRENT_LIST_DIR}/../fixtures/synthetic/cif/fixed_ions.cif")
