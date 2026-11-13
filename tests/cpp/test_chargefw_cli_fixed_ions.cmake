include("${CMAKE_CURRENT_LIST_DIR}/support/test_cli.cmake")

set(input_path "${CHARGEFW_INPUT}")
run_cli(fixed_ion_applicability 0 applicability --method sqeqp
        --parameter-set SQEqp_Schindler2021_CCD_gen
        --fixed-ions MG --fixed-ions CA "${input_path}")

fresh_output_directory(output_directory fixed_ions)
set(output_json "${output_directory}/fixed_ions.chargefw.json")
run_cli(fixed_ion_calculation 0 calculate --method sqeqp
        --parameter-set SQEqp_Schindler2021_CCD_gen --execution full
        --fixed-ions MG --fixed-ions MG --fixed-ions CA
        "${input_path}" "${output_directory}")
expect_files_exist(fixed_ion_calculation "${output_json}")
file(READ "${output_json}" result_json)
expect_json(fixed_ion_calculation "${result_json}" success results 0 status)

string(JSON charges GET "${result_json}" results 0 assignments 0 charges)
string(JSON charge_count LENGTH "${charges}")
string(JSON mg_charge GET "${charges}" 3)
string(JSON ca_charge GET "${charges}" 4)
string(JSON total_charge GET "${result_json}" results 0 assignments 0 total_charge)
if(NOT charge_count EQUAL 5 OR NOT mg_charge EQUAL 2 OR NOT ca_charge EQUAL 2 OR
   NOT total_charge EQUAL 4)
    message(FATAL_ERROR "fixed source charges or the active total are wrong: ${charges}")
endif()

string(JSON fixed_groups GET "${result_json}" calculation_provenance effective fixed_ions)
string(JSON group_count LENGTH "${fixed_groups}")
if(NOT group_count EQUAL 2)
    message(FATAL_ERROR "missing effective fixed-charge provenance: ${fixed_groups}")
endif()
expect_json(fixed_ion_provenance "${fixed_groups}" MG 0 component_id)
expect_json(fixed_ion_provenance "${fixed_groups}" CA 1 component_id)
string(JSON mg_index GET "${fixed_groups}" 0 instances 0 atom_index)
string(JSON mg_value GET "${fixed_groups}" 0 charge)
string(JSON ca_index GET "${fixed_groups}" 1 instances 0 atom_index)
string(JSON ca_value GET "${fixed_groups}" 1 charge)
if(NOT mg_index EQUAL 3 OR NOT mg_value EQUAL 2 OR NOT ca_index EQUAL 4 OR NOT ca_value EQUAL 2)
    message(FATAL_ERROR "fixed sources differ from input order or preset values: ${fixed_groups}")
endif()

run_cli(unknown_fixed_ion 2 applicability --fixed-ions UNKNOWN "${input_path}")
run_cli(missing_fixed_ion 2 applicability "${input_path}" --fixed-ions)
run_cli(several_fixed_ions_in_one_occurrence 2 applicability --fixed-ions MG CA "${input_path}")

file(REMOVE_RECURSE "${output_directory}")
