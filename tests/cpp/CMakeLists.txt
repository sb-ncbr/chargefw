add_library(chargefw_test_support
        support/test_molecules.cpp
        support/test_molecules.h
)

target_link_libraries(chargefw_test_support
        PUBLIC
        chargefw::core
)

target_include_directories(chargefw_test_support
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

if(CHARGEFW_ENABLE_WARNINGS)
    chargefw_enable_warnings(chargefw_test_support)
    # Tests use partial designated initialization to emphasize only fields relevant to each fixture.
    target_compile_options(chargefw_test_support PRIVATE -Wno-missing-field-initializers)
endif()

function(chargefw_add_snitch_test name)
    add_executable(${name}
            ${ARGN}
    )

    target_link_libraries(${name}
            PRIVATE
            chargefw::core
            chargefw_test_support
            snitch::snitch
    )

    if(CHARGEFW_ENABLE_CLANG_TIDY)
        # Snitch's expression decomposition intentionally forms an intermediate chained comparison.
        # Suppress the resulting false positive only for Snitch test targets; project targets retain
        # this clang-tidy check.
        set_property(TARGET ${name} APPEND PROPERTY CXX_CLANG_TIDY
                "--checks=-bugprone-chained-comparison")
    endif()

    target_link_options(${name} PRIVATE "LINKER:--disable-new-dtags")

    target_compile_definitions(${name}
            PRIVATE
            CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
            CHARGEFW_TEST_PARAMETER_DIR="${PROJECT_SOURCE_DIR}/data/parameters"
    )

    if(CHARGEFW_ENABLE_WARNINGS)
        chargefw_enable_warnings(${name})
        # Snitch 1.3.2 uses __COUNTER__ in its public macros. Clang 22 diagnoses each expansion at
        # the project call site under -Wpedantic even though the macro definition is a system header.
        target_compile_options(${name} PRIVATE
                -Wno-c2y-extensions
                -Wno-missing-field-initializers
        )
    endif()

    add_test(
            NAME ${name}
            COMMAND ${name}
    )
endfunction()

if(CHARGEFW_BUILD_CLI)
    set(chargefw_cli_install_prefix "${PROJECT_BINARY_DIR}/tests/cpp/chargefw_cli_prefix")
    set(chargefw_cli "${chargefw_cli_install_prefix}/${CMAKE_INSTALL_BINDIR}/chargefw")

    add_test(
            NAME test_chargefw_cli_install
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_BUILD_DIR=${PROJECT_BINARY_DIR}
                    -DCHARGEFW_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                    -DCHARGEFW_INSTALL_PREFIX=${chargefw_cli_install_prefix}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/prepare_moved_install.cmake
    )
    set_tests_properties(test_chargefw_cli_install PROPERTIES
            FIXTURES_SETUP chargefw_cli_install
    )

    add_test(
            NAME test_chargefw_cli_arguments
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_CLI=${chargefw_cli}
                    -DCHARGEFW_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/water.sdf
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_arguments.cmake
    )
    add_test(
            NAME test_chargefw_cli_outputs
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_CLI=${chargefw_cli}
                    -DCHARGEFW_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/water.sdf
                    -DCHARGEFW_MOL2_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/mol2/aromatic.mol2
                    -DCHARGEFW_MIXED_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/mixed_v2000_v3000.sdf
                    -DCHARGEFW_MALFORMED_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/malformed_then_water.sdf
                    -DCHARGEFW_VERSION=${CHARGEFW_VERSION}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_outputs.cmake
    )
    add_test(
            NAME test_chargefw_cli_structural_outputs
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_CLI=${chargefw_cli}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_structural_outputs.cmake
    )
    add_test(
            NAME test_chargefw_cli_fixed_ions
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_CLI=${chargefw_cli}
                    -DCHARGEFW_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/cif/fixed_ions.cif
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_fixed_ions.cmake
    )
    set_tests_properties(
            test_chargefw_cli_arguments
            test_chargefw_cli_outputs
            test_chargefw_cli_structural_outputs
            test_chargefw_cli_fixed_ions
            PROPERTIES FIXTURES_REQUIRED chargefw_cli_install
    )
endif()

# Keep nested builds aligned with the parent without copying its output paths or build targets.
set(chargefw_install_test_variables
        CMAKE_BUILD_TYPE
        CMAKE_C_COMPILER
        CMAKE_C_COMPILER_LAUNCHER
        CMAKE_CXX_COMPILER
        CMAKE_CXX_COMPILER_LAUNCHER
        CHARGEFW_USE_SYSTEM_DEPENDENCIES
        CHARGEFW_ENABLE_CCACHE
        CHARGEFW_ENABLE_WARNINGS
        CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS
        CHARGEFW_ENABLE_ASAN
        CHARGEFW_ENABLE_UBSAN
        CHARGEFW_ENABLE_CLANG_TIDY
)
set(chargefw_install_test_configurations Debug Release RelWithDebInfo MinSizeRel
        ${CMAKE_BUILD_TYPE})
list(REMOVE_DUPLICATES chargefw_install_test_configurations)
foreach(variable CMAKE_C_FLAGS CMAKE_CXX_FLAGS CMAKE_EXE_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS
        CMAKE_MODULE_LINKER_FLAGS CMAKE_STATIC_LINKER_FLAGS CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    list(APPEND chargefw_install_test_variables ${variable})
    foreach(configuration IN LISTS chargefw_install_test_configurations)
        string(TOUPPER "${configuration}" configuration_upper)
        list(APPEND chargefw_install_test_variables "${variable}_${configuration_upper}")
    endforeach()
endforeach()
set(chargefw_install_test_settings "")
foreach(variable IN LISTS chargefw_install_test_variables)
    if(DEFINED ${variable})
        string(APPEND chargefw_install_test_settings
                "set(${variable} [==[${${variable}}]==] CACHE STRING \"Parent build setting\" FORCE)\n")
    endif()
endforeach()
set(chargefw_install_test_cache "${CMAKE_CURRENT_BINARY_DIR}/parent_build_settings.cmake")
file(CONFIGURE OUTPUT "${chargefw_install_test_cache}"
        CONTENT "${chargefw_install_test_settings}" @ONLY)

add_test(
        NAME test_chargefw_cmake_consumer
        COMMAND ${CMAKE_COMMAND}
                -DCHARGEFW_BUILD_DIR=${PROJECT_BINARY_DIR}
                -DCHARGEFW_PARENT_CACHE=${chargefw_install_test_cache}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cmake_consumer.cmake
)
add_test(
        NAME test_chargefw_custom_install_layout
        COMMAND ${CMAKE_COMMAND}
                -DCHARGEFW_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                -DCHARGEFW_BINARY_DIR=${PROJECT_BINARY_DIR}
                -DCHARGEFW_PARENT_CACHE=${chargefw_install_test_cache}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_custom_install_layout.cmake
)

add_subdirectory(core)
add_subdirectory(charges)
add_subdirectory(calculation)
add_subdirectory(adapters)
add_subdirectory(features)
add_subdirectory(methods)
add_subdirectory(parameters)
