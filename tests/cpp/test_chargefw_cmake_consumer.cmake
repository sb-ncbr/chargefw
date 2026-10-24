set(CHARGEFW_INSTALL_PREFIX
    "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_prefix")
if(NOT DEFINED CHARGEFW_PARENT_CACHE)
    message(FATAL_ERROR "Downstream consumer test requires the parent settings path")
endif()
include("${CMAKE_CURRENT_LIST_DIR}/prepare_moved_install.cmake")

# Builds from downloaded dependencies install the public dependency packages with ChargeFW. Builds
# from installed packages rely on the packages that the ChargeFW build found, from the parent settings.
file(GLOB_RECURSE nlohmann_configurations
     "${CHARGEFW_INSTALL_PREFIX}/*nlohmann_json*Config.cmake"
     "${CHARGEFW_INSTALL_PREFIX}/*nlohmann_json*-config.cmake")
if(CHARGEFW_FETCH_DEPENDENCIES AND NOT nlohmann_configurations)
    message(FATAL_ERROR "The installed SDK does not contain the nlohmann/json CMake package")
endif()
set(nlohmann_json_arguments)
if(nlohmann_configurations)
    list(GET nlohmann_configurations 0 nlohmann_configuration)
    get_filename_component(nlohmann_json_directory "${nlohmann_configuration}" DIRECTORY)
    list(APPEND nlohmann_json_arguments -Dnlohmann_json_DIR=${nlohmann_json_directory})
endif()

file(GLOB_RECURSE gemmi_configurations
     "${CHARGEFW_INSTALL_PREFIX}/*gemmi*Config.cmake"
     "${CHARGEFW_INSTALL_PREFIX}/*gemmi*-config.cmake")
set(gemmi_arguments)
if(gemmi_configurations)
    list(GET gemmi_configurations 0 gemmi_configuration)
    get_filename_component(gemmi_directory "${gemmi_configuration}" DIRECTORY)
    list(APPEND gemmi_arguments -Dgemmi_DIR=${gemmi_directory})
endif()

set(consumer_build_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_build")
file(REMOVE_RECURSE "${consumer_build_directory}")

run_checked("Downstream ChargeFW consumer configuration"
            "${CMAKE_COMMAND}"
            -S "${CMAKE_CURRENT_LIST_DIR}/consumer"
            -B "${consumer_build_directory}"
            -C "${CHARGEFW_PARENT_CACHE}"
            -DCMAKE_PREFIX_PATH=${CHARGEFW_INSTALL_PREFIX}
            ${nlohmann_json_arguments}
            ${gemmi_arguments})
run_checked("Downstream ChargeFW consumer build"
            "${CMAKE_COMMAND}" --build "${consumer_build_directory}")
run_checked("Downstream ChargeFW consumer" "${consumer_build_directory}/chargefw_consumer")
