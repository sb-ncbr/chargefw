chargefw_add_snitch_test(test_methods
        test_builtin_methods.cpp
        test_calculation_input.cpp
        test_charge2.cpp
        test_delre.cpp
        test_electronegativity_equalization.cpp
        test_gdac.cpp
        test_method_applicability.cpp
        test_method_calculation.cpp
        test_method_options.cpp
        test_method_prerequisites.cpp
        test_method_registry.cpp
        test_mpeoe.cpp
        test_peoe.cpp
        test_sqe.cpp
        test_sqeqp.cpp
        test_tsef.cpp
        test_veem.cpp
)
target_include_directories(test_methods PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_dense_solve
        test_dense_solve.cpp
)
target_include_directories(test_dense_solve PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(test_dense_solve PRIVATE
        Eigen3::Eigen
        TBB::tbb
)
# Passes Eigen matrices to the library, so their alignment must match: it follows the instruction
# set, which these options select.
chargefw_enable_optimizations(test_dense_solve)
