chargefw_add_snitch_test(test_calculation
        test_calculation.cpp
)

chargefw_add_snitch_test(test_fixed_ion_planning
        test_fixed_ion_planning.cpp
)

chargefw_add_snitch_test(test_fixed_charge_partition
        test_fixed_charge_partition.cpp
        test_full_fixed_ions_execution.cpp
)
target_include_directories(test_fixed_charge_partition PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_observer
        test_observer.cpp
)
target_include_directories(test_observer PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_execution_policy
        test_execution_policy.cpp
)

chargefw_add_snitch_test(test_calculation_targets
        test_calculation_targets.cpp
)
target_include_directories(test_calculation_targets PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(test_calculation_targets PRIVATE
        TBB::tbb
)

chargefw_add_snitch_test(test_reduced_execution
        test_reduced_execution.cpp
)
target_include_directories(test_reduced_execution PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_representative_execution
        test_representative_execution.cpp
)

chargefw_add_snitch_test(test_cover_execution
        test_cover_execution.cpp
)
target_include_directories(test_cover_execution PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_derivatives
        test_derivatives.cpp
)
