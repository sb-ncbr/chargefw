# The command-line sources are compiled once. The native executable links the shared library, and the
# Python package links them into its own copy of the command line (python/CMakeLists.txt).
add_library(chargefw_cli_objects OBJECT
        main.cpp
        input.cpp
        reports.cpp
        output.cpp
)

chargefw_enable_optimizations(chargefw_cli_objects)

if(CHARGEFW_ENABLE_WARNINGS)
    chargefw_enable_warnings(chargefw_cli_objects)
    # The CLI deliberately uses partial designated initialization for public aggregate result types;
    # omitted fields retain their default member initializers.
    if(NOT MSVC)
        target_compile_options(chargefw_cli_objects PRIVATE -Wno-missing-field-initializers)
    endif()
endif()

target_link_libraries(chargefw_cli_objects
        PRIVATE
        $<COMPILE_ONLY:chargefw_core_objects>
        chargefw_platform
        CLI11::CLI11
)

target_include_directories(chargefw_cli_objects PRIVATE "${PROJECT_SOURCE_DIR}/src")

# Both executables run with UTF-8 as their Windows code page; MSVC embeds a manifest listed among an
# executable's sources.
if(MSVC)
    target_sources(chargefw_cli_objects INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/utf8.manifest")
endif()

if(CHARGEFW_BUILD_CLI)
    add_executable(chargefw_cli)

    set_target_properties(chargefw_cli
            PROPERTIES
            OUTPUT_NAME chargefw
    )

    target_link_libraries(chargefw_cli
            PRIVATE
            chargefw_cli_objects
            chargefw::core
    )

    install(TARGETS chargefw_cli
            RUNTIME DESTINATION bin
    )
endif()
