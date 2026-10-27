add_subdirectory(io)
add_subdirectory(models)
add_subdirectory(classification)

target_sources(chargefw_core_objects
        PRIVATE
        parameter_helpers.h
        parameter_helpers.cpp
)


