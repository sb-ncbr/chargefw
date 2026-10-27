# Operating-system services used by the library and the command line. Each platform directory
# implements the headers in this directory under the same file names; this is the only place that
# selects a platform.
add_library(chargefw_platform STATIC)
# Linked into the shared library, but not part of its interface.
set_target_properties(chargefw_platform
        PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
)

target_compile_features(chargefw_platform PUBLIC cxx_std_23)
target_include_directories(chargefw_platform PUBLIC ${PROJECT_SOURCE_DIR}/src)

if(WIN32)
    set(platform_directory windows)
    target_link_libraries(chargefw_platform PRIVATE psapi)
else()
    set(platform_directory unix)
    target_link_libraries(chargefw_platform PRIVATE ${CMAKE_DL_LIBS})
endif()

target_sources(chargefw_platform
        PRIVATE
        ${platform_directory}/library_location.cpp
        ${platform_directory}/peak_memory.cpp
)

if(CHARGEFW_ENABLE_WARNINGS)
    chargefw_enable_warnings(chargefw_platform)
endif()

if(CHARGEFW_ENABLE_ASAN)
    chargefw_enable_sanitizer(chargefw_platform address)
endif()

if(CHARGEFW_ENABLE_UBSAN)
    chargefw_enable_sanitizer(chargefw_platform undefined)
endif()
