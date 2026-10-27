cmake_minimum_required(VERSION 3.28)

project(
  chargefw
  VERSION 3.0.0
  DESCRIPTION "Empirical partial atomic charge calculation framework"
  LANGUAGES CXX
)

# PEP 440 pre-release suffix (aN, bN, or rcN); empty for final releases.
set(CHARGEFW_VERSION_SUFFIX "a3")
if(NOT CHARGEFW_VERSION_SUFFIX MATCHES "^((a|b|rc)(0|[1-9][0-9]*))?$")
  message(FATAL_ERROR "Invalid CHARGEFW_VERSION_SUFFIX: ${CHARGEFW_VERSION_SUFFIX}")
endif()
set(CHARGEFW_VERSION "${PROJECT_VERSION}${CHARGEFW_VERSION_SUFFIX}")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# The bundled parameter sets relative to the file that contains the library code: the shared library in
# a native installation, or the extension and command line in the Python package.
file(RELATIVE_PATH CHARGEFW_PARAMETER_DIRECTORY_FROM_LIBRARY
  "${CMAKE_INSTALL_FULL_LIBDIR}"
  "${CMAKE_INSTALL_FULL_DATADIR}/chargefw/parameters"
)
set(CHARGEFW_PARAMETER_DIRECTORY_IN_PACKAGE "_data/parameters")

option(CHARGEFW_ENABLE_CCACHE "Use ccache if available" ON)
option(CHARGEFW_BUILD_TESTS "Build tests" ${PROJECT_IS_TOP_LEVEL})
option(CHARGEFW_BUILD_CLI "Build the chargefw command-line application" ${PROJECT_IS_TOP_LEVEL})
option(CHARGEFW_BUILD_PYTHON "Build the ChargeFW Python extension" OFF)
set(CHARGEFW_PYTHON_EXECUTABLE "" CACHE FILEPATH
    "Python interpreter for ChargeFW bindings (empty uses CMake discovery)")
option(CHARGEFW_FETCH_DEPENDENCIES
       "Download pinned dependencies instead of using installed packages" OFF)
option(CHARGEFW_ENABLE_WARNINGS "Enable project warnings" ON)
option(CHARGEFW_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(CHARGEFW_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(CHARGEFW_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
option(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS "Optimize for the build host CPU" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# The Visual Studio generator compiles the sources of one target sequentially unless /MP is set. It
# applies to the fetched dependencies as well, which are added below.
add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)

include(cmake/dependencies.cmake)
include(cmake/developer_tools.cmake)
include(cmake/optimization.cmake)
include(cmake/runtime_paths.cmake)
include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)

chargefw_setup_dependencies()
chargefw_setup_developer_tools()

# The library sources are compiled once. The installed shared library is built from these objects, and
# the Python extension links them directly, so wheels contain no separate ChargeFW library.
add_library(chargefw_core_objects OBJECT)
set_target_properties(chargefw_core_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_compile_features(chargefw_core_objects
  PUBLIC
    cxx_std_23
)

target_include_directories(chargefw_core_objects
  PRIVATE
    src

  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
)

target_link_libraries(chargefw_core_objects
  PUBLIC
    nlohmann_json::nlohmann_json
    gemmi::gemmi_cpp
  PRIVATE
    $<BUILD_LOCAL_INTERFACE:Eigen3::Eigen>
    $<BUILD_LOCAL_INTERFACE:nanoflann::nanoflann>
    TBB::tbb
    $<BUILD_LOCAL_INTERFACE:chargefw_platform>
)

target_compile_definitions(chargefw_core_objects
  PRIVATE
    EIGEN_DONT_PARALLELIZE
)

chargefw_enable_optimizations(chargefw_core_objects)

configure_file(
  ${PROJECT_SOURCE_DIR}/include/chargefw/config.h.in
  ${PROJECT_BINARY_DIR}/generated/chargefw/config.h
  @ONLY
)

if(CHARGEFW_ENABLE_WARNINGS)
  chargefw_enable_warnings(chargefw_core_objects)
endif()

add_library(chargefw_core SHARED)
add_library(chargefw::core ALIAS chargefw_core)
set_target_properties(chargefw_core PROPERTIES EXPORT_NAME core)

target_compile_features(chargefw_core
  PUBLIC
    cxx_std_23
)

target_include_directories(chargefw_core
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(chargefw_core
  PUBLIC
    nlohmann_json::nlohmann_json
    gemmi::gemmi_cpp
  PRIVATE
    $<BUILD_LOCAL_INTERFACE:chargefw_core_objects>
)

# Sanitizer options are public: the objects need them to compile, and programs linking the library need
# them to link.
foreach(target IN ITEMS chargefw_core_objects chargefw_core)
  if(CHARGEFW_ENABLE_ASAN)
    chargefw_enable_sanitizer(${target} address)
  endif()

  if(CHARGEFW_ENABLE_UBSAN)
    chargefw_enable_sanitizer(${target} undefined)
  endif()
endforeach()

add_subdirectory(src)

# The Python package contains its own copy of the command line, so its sources are needed for either.
if(CHARGEFW_BUILD_CLI OR CHARGEFW_BUILD_PYTHON)
  add_subdirectory(apps/chargefw)
endif()

if(CHARGEFW_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(CHARGEFW_BUILD_TESTS)
  chargefw_setup_test_dependencies()
  enable_testing()
  add_subdirectory(tests)
endif()

include(cmake/install_rules.cmake)
