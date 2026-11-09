chargefw_add_snitch_test(test_normal_equations
        test_normal_equations.cpp
)

chargefw_add_snitch_test(test_method_derivatives
        test_method_derivatives.cpp
)
target_include_directories(test_method_derivatives PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_parameter_invariances
        test_parameter_invariances.cpp
)

chargefw_add_snitch_test(test_parameter_template
        test_parameter_template.cpp
)

chargefw_add_snitch_test(test_parameter_vector
        test_parameter_vector.cpp
)

chargefw_add_snitch_test(test_parameter_entries
        test_parameter_entries.cpp
)
