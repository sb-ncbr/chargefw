chargefw_add_snitch_test(test_mol
        test_mol.cpp
)

chargefw_add_snitch_test(test_json
        test_json.cpp
)

chargefw_add_snitch_test(test_charge_result
        test_charge_result.cpp
)

chargefw_add_snitch_test(test_json_output
        test_json_output.cpp
)

chargefw_add_snitch_test(test_mol2_output
        test_mol2_output.cpp
)

chargefw_add_snitch_test(test_pqr_output
        test_pqr_output.cpp
)

chargefw_add_snitch_test(test_pdb
        test_pdb.cpp
)

chargefw_add_snitch_test(test_mmcif
        test_mmcif.cpp
)

chargefw_add_snitch_test(test_mmcif_output
        test_mmcif_output.cpp
)

chargefw_add_snitch_test(test_fixed_ions
        test_fixed_ions.cpp
)

target_include_directories(test_fixed_ions PRIVATE ${PROJECT_SOURCE_DIR}/src)
