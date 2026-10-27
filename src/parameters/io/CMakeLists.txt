target_sources(chargefw_core_objects
        PRIVATE
        parameter_set_io.cpp
        parameter_data_paths.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameters/io/parameter_set_io.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameters/io/parameter_data_paths.h
)