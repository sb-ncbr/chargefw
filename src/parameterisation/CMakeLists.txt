target_sources(chargefw_core_objects
        PRIVATE
        normal_equations.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameterisation/normal_equations.h
)
