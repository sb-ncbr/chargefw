target_sources(chargefw_core_objects
        PRIVATE
        eem_derivatives.cpp
        method_derivatives.cpp
        method_derivatives.h
        normal_equations.cpp
        parameter_columns.h
        sqe_derivatives.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameterisation/normal_equations.h
)
