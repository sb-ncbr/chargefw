include("${CMAKE_CURRENT_LIST_DIR}/support/test_cli.cmake")

function(run_structural_input extension contents input_stem expected_selection expected_bonds)
    set(input_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${input_stem}.${extension}")
    set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${input_stem}_outputs")
    set(output_prefix "${output_directory}/chargefw_cli_${input_stem}.chargefw")
    file(WRITE "${input_path}" "${contents}")
    file(REMOVE_RECURSE "${output_directory}")

    run_cli("${input_stem}" 0 calculate --output-mol2 --output-mmcif --output-pqr ${ARGN}
            "${input_path}" "${output_directory}")

    expect_files_exist("${input_stem}" "${output_prefix}.json" "${output_prefix}.cif"
                       "${output_prefix}.mol2" "${output_prefix}.pqr")
    file(READ "${output_prefix}.json" json_output)
    expect_json("${input_stem}" "${json_output}" success results 0 status)
    expect_json("${input_stem}" "${json_output}" "${expected_selection}"
                results 0 input import policy selection)
    expect_json("${input_stem}" "${json_output}" "${expected_bonds}"
                results 0 input import policy bonds)

    file(REMOVE "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

function(expect_non_structural_options_rejected)
    set(input_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_nonstructural.json")
    set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_nonstructural_outputs")
    file(WRITE "${input_path}" "{\"schema_version\": \"1.0\", \"molecules\": [{\"atoms\": [{\"atomic_number\": 8, \"formal_charge\": 0}]}]}\n")
    file(REMOVE_RECURSE "${output_directory}")

    run_cli(non_structural_options 2 calculate --bonds templates
            "${input_path}" "${output_directory}")

    file(REMOVE "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

# --components names the dictionary that bonds a residue without a built-in template; the result
# records the dictionary's version.
function(expect_component_dictionary)
    set(input_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_components.pdb")
    set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_components_outputs")
    file(WRITE "${input_path}" "HETATM    1  C   FOR A   1       0.000   0.000   0.000  1.00 20.00           C\nHETATM    2  O   FOR A   1       1.200   0.000   0.000  1.00 20.00           O\nHETATM    3  H1  FOR A   1      -0.500   0.900   0.000  1.00 20.00           H\nHETATM    4  H2  FOR A   1      -0.500  -0.900   0.000  1.00 20.00           H\nEND\n")
    file(REMOVE_RECURSE "${output_directory}")

    run_cli(components 0 calculate --method formal --components "${CHARGEFW_COMPONENTS}"
            "${input_path}" "${output_directory}")
    file(READ "${output_directory}/chargefw_cli_components.chargefw.json" json_output)
    expect_json(components "${json_output}" "2026.11.13.post1"
                results 0 input import policy components)

    # A file that is not a dictionary, here the input itself, is refused before reading.
    file(REMOVE_RECURSE "${output_directory}")
    run_cli(components_other_file 2 calculate --method formal --components "${input_path}"
            "${input_path}" "${output_directory}")
    if(EXISTS "${output_directory}")
        message(FATAL_ERROR "components_other_file: output was written")
    endif()

    file(REMOVE "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

# --annotate-mmcif writes the mmCIF input with charges; charges already there need --if-exists.
function(expect_annotated_mmcif contents)
    set(input_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_annotate.cif")
    set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_annotate_outputs")
    set(annotated "${output_directory}/chargefw_cli_annotate.chargefw.annotated.cif")
    file(WRITE "${input_path}" "${contents}")
    file(REMOVE_RECURSE "${output_directory}")

    run_cli(annotate 0 calculate --method formal --annotate-mmcif
            "${input_path}" "${output_directory}")
    expect_files_exist(annotate "${annotated}")
    file(READ "${annotated}" annotated_contents)
    if(NOT annotated_contents MATCHES "_note\\.text" OR
       NOT annotated_contents MATCHES "_sb_ncbr_partial_atomic_charges\\.charge")
        message(FATAL_ERROR "annotate: the annotated mmCIF lacks input data or charges")
    endif()

    # The annotated file as input: refused before calculating, so nothing is written.
    file(COPY_FILE "${annotated}" "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
    run_cli(annotate_existing 2 calculate --method formal --annotate-mmcif
            "${input_path}" "${output_directory}")
    expect_files_absent(annotate_existing "${output_directory}/chargefw_cli_annotate.chargefw.json")
    run_cli(annotate_replace 0 calculate --method formal --annotate-mmcif --if-exists replace
            "${input_path}" "${output_directory}")
    expect_files_exist(annotate_replace "${annotated}")
    run_cli(annotate_append 0 calculate --method dummy --annotate-mmcif --if-exists append
            "${input_path}" "${output_directory}")
    file(READ "${annotated}" annotated_contents)
    if(NOT annotated_contents MATCHES "formal" OR NOT annotated_contents MATCHES "dummy")
        message(FATAL_ERROR "annotate_append: the annotated mmCIF lacks one of the charge sets")
    endif()
    run_cli(annotate_option_alone 2 calculate --method formal --if-exists replace
            "${input_path}" "${output_directory}")

    # Charge metadata without charge rows is refused before calculating as well.
    file(WRITE "${input_path}" "${contents}_sb_ncbr_partial_atomic_charges_meta.id 1\n")
    file(REMOVE_RECURSE "${output_directory}")
    run_cli(annotate_metadata_only 2 calculate --method formal --annotate-mmcif
            "${input_path}" "${output_directory}")
    expect_files_absent(annotate_metadata_only
                        "${output_directory}/chargefw_cli_annotate.chargefw.json")

    # A PDB input is converted; an input without a structure cannot be annotated.
    set(pdb_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_annotate.pdb")
    file(WRITE "${pdb_path}" "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  \nEND\n")
    file(REMOVE_RECURSE "${output_directory}")
    run_cli(annotate_pdb 0 calculate --method formal --annotate-mmcif
            "${pdb_path}" "${output_directory}")
    expect_files_exist(annotate_pdb "${annotated}")
    set(json_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_annotate.json")
    file(WRITE "${json_path}" "{\"schema_version\": \"1.0\", \"molecules\": [{\"atoms\": [{\"atomic_number\": 8, \"formal_charge\": 0}]}]}\n")
    run_cli(annotate_json 2 calculate --method formal --annotate-mmcif
            "${json_path}" "${output_directory}")
    file(REMOVE "${json_path}")

    file(REMOVE "${input_path}" "${pdb_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

# A .gz file is read as the format of the extension before .gz, and outputs are named without both.
function(expect_gzip_input contents)
    set(plain_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_gzip.cif")
    set(input_path "${plain_path}.gz")
    fresh_output_directory(output_directory gzip_outputs)
    file(WRITE "${plain_path}" "${contents}")
    file(ARCHIVE_CREATE OUTPUT "${input_path}" PATHS "${plain_path}" FORMAT raw COMPRESSION GZip)

    run_cli(gzip 0 calculate --method formal --output-mmcif "${input_path}" "${output_directory}")
    expect_files_exist(gzip "${output_directory}/chargefw_cli_gzip.chargefw.json"
                       "${output_directory}/chargefw_cli_gzip.chargefw.cif")

    # Text named .gz is not gzip data.
    file(COPY_FILE "${plain_path}" "${input_path}")
    fresh_output_directory(output_directory gzip_outputs)
    run_cli(gzip_uncompressed 2 calculate --method formal "${input_path}" "${output_directory}")
    expect_files_absent(gzip_uncompressed "${output_directory}/chargefw_cli_gzip.chargefw.json")

    file(REMOVE "${plain_path}" "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

set(water_atom_site "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\nHETATM 1 O O . HOH A 1 ? 0.000 0.000 0.000 1.00 20.00 0 1 HOH A O 1\nHETATM 2 H H1 . HOH A 1 ? 0.957 0.000 0.000 1.00 20.00 0 1 HOH A H1 1\nHETATM 3 H H2 . HOH A 1 ? -0.239 0.927 0.000 1.00 20.00 0 1 HOH A H2 1\n#\n")

run_structural_input(
        pdb
        "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  \nATOM      2  H1  HOH A   1       0.957   0.000   0.000  1.00 20.00           H  \nATOM      3  H2  HOH A   1      -0.239   0.927   0.000  1.00 20.00           H  \nEND\n"
        structural_pdb
        all
        templates
        --selection all
        --bonds templates
)

run_structural_input(
        cif
        "data_structural_default_bonds\n${water_atom_site}"
        structural_default_bonds
        all
        hybrid
)

expect_annotated_mmcif("data_annotate\n_note.text kept\n${water_atom_site}")

expect_non_structural_options_rejected()
expect_component_dictionary()
expect_gzip_input("data_gzip\n${water_atom_site}")
