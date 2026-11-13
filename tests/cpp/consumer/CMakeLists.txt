cmake_minimum_required(VERSION 3.28)

project(chargefw_consumer LANGUAGES CXX)

find_package(chargefw CONFIG REQUIRED)

add_executable(chargefw_consumer main.cpp)
target_link_libraries(chargefw_consumer PRIVATE chargefw::core)
# The installed version header must describe the release the package reports.
target_compile_definitions(chargefw_consumer PRIVATE
    EXPECTED_VERSION_MAJOR=${chargefw_VERSION_MAJOR}
    EXPECTED_VERSION_MINOR=${chargefw_VERSION_MINOR}
    EXPECTED_VERSION_PATCH=${chargefw_VERSION_PATCH}
    EXPECTED_VERSION="${chargefw_VERSION}")
