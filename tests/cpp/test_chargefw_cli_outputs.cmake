include("${CMAKE_CURRENT_LIST_DIR}/support/test_cli.cmake")

fresh_output_directory(json_failure_directory json_failure)
file(MAKE_DIRECTORY "${json_failure_directory}/water.chargefw.json")
run_cli(json_output_failure 2 calculate --method eem "${CHARGEFW_INPUT}"
        "${json_failure_directory}")
file(REMOVE_RECURSE "${json_failure_directory}")

fresh_output_directory(export_failure_directory export_failure)
set(export_failure_prefix "${export_failure_directory}/water.chargefw")
file(MAKE_DIRECTORY "${export_failure_prefix}.mol2")
run_cli(molecular_export_failure 6 calculate --method eem --output-mol2 "${CHARGEFW_INPUT}"
        "${export_failure_directory}")
expect_files_exist(molecular_export_failure "${export_failure_prefix}.json")
file(READ "${export_failure_prefix}.json" export_failure_json)
expect_json(molecular_export_failure "${export_failure_json}" success status)
file(REMOVE_RECURSE "${export_failure_directory}")

foreach(mode IN ITEMS full cutoff cover)
    fresh_output_directory(mode_directory "${mode}")
    set(execution_arguments --method eem --execution ${mode})
    if(NOT mode STREQUAL "full")
        list(APPEND execution_arguments --radius 8)
    endif()
    if(mode STREQUAL "cutoff")
        list(APPEND execution_arguments --progress)
    endif()

    run_cli("${mode}" 0 calculate ${execution_arguments} "${CHARGEFW_INPUT}" "${mode_directory}")
    file(READ "${mode_directory}/water.chargefw.json" mode_json)
    expect_json("${mode}" "${mode_json}" success results 0 status)
    expect_json("${mode}" "${mode_json}" "${mode}"
                calculation_provenance effective execution mode)
    file(REMOVE_RECURSE "${mode_directory}")
endforeach()

function(expect_invalid_policy label)
    fresh_output_directory(policy_directory "invalid_${label}")
    run_cli("${label}" 2 calculate ${ARGN} "${CHARGEFW_INPUT}" "${policy_directory}")
    expect_files_absent("${label}" "${policy_directory}/water.chargefw.json")
    file(REMOVE_RECURSE "${policy_directory}")
endfunction()

expect_invalid_policy(missing_radius --execution cutoff)
expect_invalid_policy(short_radius --execution cover --radius 7)
expect_invalid_policy(parameter_without_method --parameter-set QEq_original)

fresh_output_directory(malformed_directory malformed)
run_cli(malformed_input 2 calculate "${CHARGEFW_MALFORMED_INPUT}" "${malformed_directory}")
expect_files_absent(malformed_input "${malformed_directory}/malformed_then_water.chargefw.json")
file(REMOVE_RECURSE "${malformed_directory}")

fresh_output_directory(mixed_directory mixed)
run_cli(mixed_records 0 calculate --method eem --execution full "${CHARGEFW_MIXED_INPUT}"
        "${mixed_directory}")
file(READ "${mixed_directory}/mixed_v2000_v3000.chargefw.json" mixed_json)
string(JSON mixed_record_count LENGTH "${mixed_json}" results)
if(NOT mixed_record_count EQUAL 2)
    message(FATAL_ERROR "mixed V2000/V3000 records were not preserved")
endif()
expect_json(mixed_records "${mixed_json}" v2000 results 0 input record_id)
expect_json(mixed_records "${mixed_json}" v3000 results 1 input record_id)
file(REMOVE_RECURSE "${mixed_directory}")

set(multiconformer_input "${CMAKE_CURRENT_BINARY_DIR}/multi_conformer.json")
file(WRITE "${multiconformer_input}" [=[
{
  "schema_version": "1.0",
  "molecules": [{
    "id": "oxygen",
    "atoms": [{"atomic_number": 8, "formal_charge": 0}],
    "conformers": [
      {"id": "first", "coordinates": [[1.25, 2.5, 3.75]]},
      {"id": "second", "coordinates": [[9.5, 8.5, 7.5]]}
    ]
  }]
}
]=])
fresh_output_directory(multiconformer_directory multiconformer)
set(multiconformer_prefix "${multiconformer_directory}/multi_conformer.chargefw")
run_cli(multiconformer 0 calculate --method formal --output-mol2 --output-mmcif
        "${multiconformer_input}" "${multiconformer_directory}")
expect_files_exist(multiconformer "${multiconformer_prefix}.json" "${multiconformer_prefix}.mol2"
                   "${multiconformer_prefix}.cif")
file(REMOVE "${multiconformer_input}")
file(REMOVE_RECURSE "${multiconformer_directory}")

set(coordinate_free_input "${CMAKE_CURRENT_BINARY_DIR}/coordinate_free.json")
file(WRITE "${coordinate_free_input}" [=[
{"schema_version":"1.0","molecules":[{"atoms":[{"atomic_number":1,"formal_charge":0}]}]}
]=])
fresh_output_directory(coordinate_free_directory coordinate_free)
file(MAKE_DIRECTORY "${coordinate_free_directory}")
set(refused_mol2 "${coordinate_free_directory}/coordinate_free.chargefw.mol2")
# An unrelated file that looks like a temporary output must survive every export.
set(unrelated_temporary "${coordinate_free_directory}/coordinate_free.chargefw.mol2.tmp")
file(WRITE "${unrelated_temporary}" "unrelated\n")
function(expect_refused_export label expected_mol2)
  run_cli("${label}" 6 calculate --method formal --output-mol2 "${coordinate_free_input}"
          "${coordinate_free_directory}")
  expect_files_exist("${label}" "${coordinate_free_directory}/coordinate_free.chargefw.json")
  if(expected_mol2 STREQUAL "")
    expect_files_absent("${label}" "${refused_mol2}")
  else()
    file(READ "${refused_mol2}" mol2_contents)
    if(NOT mol2_contents STREQUAL expected_mol2)
      message(FATAL_ERROR "${label}: refused MOL2 export changed the existing file")
    endif()
  endif()
  file(READ "${unrelated_temporary}" unrelated_contents)
  file(GLOB temporaries "${coordinate_free_directory}/*.tmp")
  if(NOT unrelated_contents STREQUAL "unrelated\n" OR NOT temporaries STREQUAL unrelated_temporary)
    message(FATAL_ERROR "${label}: temporary files left or changed: ${temporaries}")
  endif()
endfunction()
# A refused export creates no file, and leaves an earlier file of the same name unchanged.
expect_refused_export(coordinate_free_export "")
file(WRITE "${refused_mol2}" "previous export\n")
expect_refused_export(coordinate_free_existing_export "previous export\n")
file(REMOVE "${coordinate_free_input}")
file(REMOVE_RECURSE "${coordinate_free_directory}")

set(range_input "${CMAKE_CURRENT_BINARY_DIR}/out_of_range_charge.json")
file(WRITE "${range_input}" [=[
{"schema_version":"1.0","molecules":[{"atoms":[{"atomic_number":1,"formal_charge":6}],"conformers":[{"coordinates":[[0,0,0]]}]}]}
]=])
fresh_output_directory(range_directory out_of_range_charge)
set(range_prefix "${range_directory}/out_of_range_charge.chargefw")
run_cli(out_of_range_export 6 calculate --method formal --output-mol2 --output-mmcif
        "${range_input}" "${range_directory}")
expect_files_exist(out_of_range_export "${range_prefix}.json" "${range_prefix}.mol2")
file(REMOVE "${range_input}")
file(REMOVE_RECURSE "${range_directory}")

fresh_output_directory(no_plan_directory no_plan)
set(no_plan_prefix "${no_plan_directory}/water.chargefw")
run_cli(no_plan 3 calculate --method smpqeq "${CHARGEFW_INPUT}" "${no_plan_directory}")
expect_files_exist(no_plan "${no_plan_prefix}.json")
expect_files_absent(no_plan "${no_plan_prefix}.mol2" "${no_plan_prefix}.cif")
file(READ "${no_plan_prefix}.json" no_plan_json)
expect_json(no_plan "${no_plan_json}" no_executable_plan status)
expect_json(no_plan "${no_plan_json}" no_executable_plan results 0 status)
expect_json(no_plan "${no_plan_json}" no_executable_plan diagnostics 0 code)
string(JSON no_plan_cause GET "${no_plan_json}" results 0 diagnostics 1 code)
if(no_plan_cause STREQUAL "no_executable_plan")
    message(FATAL_ERROR "no-plan CLI JSON does not report the rejection cause")
endif()
file(REMOVE_RECURSE "${no_plan_directory}")

fresh_output_directory(import_warning_directory input_warning)
run_cli(import_warning 0 calculate --method eem "${CHARGEFW_MOL2_INPUT}"
        "${import_warning_directory}")
file(READ "${import_warning_directory}/aromatic.chargefw.json" import_warning_json)
expect_json(import_warning "${import_warning_json}" partial_charges_ignored
            results 0 diagnostics 0 code)
file(REMOVE_RECURSE "${import_warning_directory}")

fresh_output_directory(invalid_directory invalid)
set(invalid_json_path "${invalid_directory}/water.chargefw.json")
run_cli(invalid_request 2 calculate --cutoff-atom-threshold unlimited --cover-atom-threshold 1
        "${CHARGEFW_INPUT}" "${invalid_directory}")
expect_files_exist(invalid_request "${invalid_json_path}")
file(READ "${invalid_json_path}" invalid_json)
expect_json(invalid_request "${invalid_json}" invalid_input_or_request status)
expect_json(invalid_request "${invalid_json}" invalid_input_or_request diagnostics 0 code)
file(REMOVE_RECURSE "${invalid_directory}")

fresh_output_directory(default_directory outputs)
set(default_prefix "${default_directory}/water.chargefw")
run_cli(default_output 0 calculate "${CHARGEFW_INPUT}" "${default_directory}")
expect_files_exist(default_output "${default_prefix}.json")
expect_files_absent(default_output "${default_prefix}.mol2" "${default_prefix}.cif")
file(READ "${default_prefix}.json" default_json)
expect_json(default_output "${default_json}" 1.0 schema_version)
expect_json(default_output "${default_json}" "${CHARGEFW_VERSION}" generator version)
expect_json(default_output "${default_json}" success results 0 status)
expect_json(default_output "${default_json}" auto calculation_provenance requested execution kind)
expect_json(default_output "${default_json}" full calculation_provenance effective execution mode)
expect_json(default_output "${default_json}" 20000
            calculation_provenance requested resource_policy cutoff_atom_threshold)
expect_json(default_output "${default_json}" 80000
            calculation_provenance requested resource_policy cover_atom_threshold)
string(JSON selected_method GET "${default_json}" calculation_provenance effective method id)
if(selected_method STREQUAL "")
    message(FATAL_ERROR "Selected method provenance is empty")
endif()
string(JSON permissive_types GET "${default_json}"
       calculation_provenance requested classification permissive_types)
if(permissive_types)
    message(FATAL_ERROR "Expected strict type classification by default")
endif()
file(REMOVE_RECURSE "${default_directory}")

fresh_output_directory(explicit_directory explicit_outputs)
set(explicit_prefix "${explicit_directory}/water.chargefw")
run_cli(explicit_outputs 0 calculate --output-mol2 --output-mmcif "${CHARGEFW_INPUT}"
        "${explicit_directory}")
expect_files_exist(explicit_outputs "${explicit_prefix}.json" "${explicit_prefix}.mol2"
                   "${explicit_prefix}.cif")
file(REMOVE_RECURSE "${explicit_directory}")
