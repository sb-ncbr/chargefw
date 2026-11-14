set(CHARGEFW_INSTALL_PREFIX
    "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_prefix")
if(NOT DEFINED CHARGEFW_PARENT_CACHE)
    message(FATAL_ERROR "Downstream consumer test requires the parent settings path")
endif()
include("${CMAKE_CURRENT_LIST_DIR}/prepare_moved_install.cmake")

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
