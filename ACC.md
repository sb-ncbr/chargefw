# Moving Atomic Charge Calculator III to ChargeFW 3

Working notes for the migration of [Atomic Charge Calculator III](https://acc.biodata.ceitec.cz/)
(ACC) from ChargeFW2 to ChargeFW 3, and for the ChargeFW changes the migration calls for. They record an
analysis made on 2026-11-14 with ChargeFW 3.0.0rc5. Nothing in ACC has been changed yet. Remove this
file when the migration is complete; move what remains open to `TODO.md`.

## Goal and rules

- ACC users should not notice the switch. The published behavior of ACC is the reference.
- Bug fixes and improvements are acceptable; new features are not a goal. The maintainer decides each
  visible difference; the section [Decisions](#decisions) lists what is decided and what is open.
- Mol* rendering of results must not get worse in any case.
- ACC work happens only on the `dev` branch of ACC.
- Where a mitigation would otherwise be needed in every consumer, fix ChargeFW instead of ACC.

## Where things are

| What | Where |
| --- | --- |
| ACC checkout | `~/Projects/AtomicChargeCalculator`, branch `dev` (`origin/master` and `origin/dev` were the same commit, `d9e537f`) |
| ChargeFW2 as ACC uses it | submodule `ChargeFW2` at `ce2baa0` (branch `dev` of `sb-ncbr/ChargeFW2`), built into the image `chargefw2-base:acc` |
| ChargeFW2 Python bindings | `ChargeFW2/src/python.cpp` |
| ACC backend, ChargeFW2 calls | `src/backend/app/integrations/chargefw2/`, `src/backend/app/services/chargefw2.py` |
| ACC backend, viewer files | `src/backend/app/services/mmcif.py` |
| ACC frontend, viewer | `src/frontend/acc/src/lib/viewer/viewer.ts`, `src/frontend/acc/src/components/results/` |
| ACC examples | `src/backend/app/examples/` (outputs of ChargeFW2, not the inputs) |

The untracked `chargefw.md` in the ACC checkout is an earlier note from 2026-10-29. It predates the
release candidates and is superseded by this file.

## How ACC uses ChargeFW2

ACC calls six functions of the `chargefw2` module from a pool of four threads:

| ChargeFW2 call | What ACC does with it |
| --- | --- |
| `Molecules(file, read_hetatm, ignore_water, permissive_types)` | Reads an uploaded SDF, MOL2, PDB, or mmCIF file. The three flags are ACC's advanced settings; defaults `True`, `False`, `True` |
| `molecules.info()` | Molecule count, atom count, and atom counts per element, shown for each file |
| `get_available_methods()`, `get_available_parameters(method)` | Method and parameter catalogs |
| `get_suitable_methods(molecules)` | Methods with their applicable parameter sets, in priority order. ACC intersects them over all files of a calculation. The first method and its first set are the default model |
| `get_best_parameters(molecules, method, permissive_types)` | The first applicable set of one method |
| `calculate_charges(molecules, method, parameters, out_dir)` | Returns `{molecule name: [charges]}` and writes the output files |

Behavior of ChargeFW2 that ACC inherits:

- `calculate_charges` writes into `out_dir`: one `<name>.fw2.cif` per molecule (the name in lower case),
  `<input>-<method>-<parameters>.<ext>.txt`, and either a `.pqr` (structures) or a `.mol2` file.
- Molecule names are sanitized: every character other than an ASCII letter or digit becomes `_`, the
  name is cut at 80 characters, and duplicates get `_0`, `_1`, … A structure is named by `_entry.id`,
  from the PDB `HEADER` or the mmCIF. ACC uses the names as dictionary keys, file names, and the
  `molecule` URL parameter; the TXT file has them in upper case.
- For an mmCIF upload, the `.fw2.cif` is the uploaded data block with the atoms that were not read
  removed, `_chem_comp` removed, and the charge categories added. For a PDB upload, Gemmi converts the
  file to mmCIF first (`setup_entities`, `assign_label_seq_id`, `make_mmcif_block`). For SDF and MOL2,
  ChargeFW2 writes `_atom_site` and `_chem_comp_bond` itself.
- ACC then replaces the charge categories in each `.fw2.cif` with its own
  (`MmCIFService._write_molecule_to_mmcif`): one row per model in
  `_sb_ncbr_partial_atomic_charges_meta` with `method` set to `'<method>/<parameters>'`, charges
  dictionary 1.0, and `atom_id` equal to the atom's position plus one.
- Only the first model of a file is read. Atoms are kept when they have no alternate location or
  alternate location `A`.
- `read_hetatm=False` drops every `HETATM` record, modified residues of the polymer included.
  `ignore_water=True` drops residues named `HOH`.
- Bonds of structures come from residue templates only: `amino_acids.txt` and `other_residues.txt`, a
  23 MB copy of the bonds of the whole Chemical Component Dictionary (CCD), matched by residue and
  atom names. Bonds in the file (`CONECT`, `_chem_comp_bond`) are ignored. Backbone bonds join `C` of
  residue number *n* to `N` of residue number *n* + 1.
- Structures above 20,000 atoms get only the methods that support a cutoff, which then run with a
  12 Å cutoff.
- A molecule whose charges are not finite is left out of the result without an error.
- The reading options are global variables of the library, set by each `Molecules(...)` call. With
  ACC's four threads this is a race; ChargeFW 3 has no such state.

ACC stores calculations in its database by file hash, method, parameter-set ID, and the three
settings, and reuses a stored result for the same request (`filter_existing_calculations`).

## How the comparison was made

Both libraries ran on the same inputs with ACC's default settings. ChargeFW2 ran in a container made
from `chargefw2-base:acc` plus `python3`, through the same bindings ACC uses. ChargeFW ran from the
3.0.0rc5 wheel from PyPI with `method_scope="all"`, `parameter_matching="permissive"`, and
`selection` translated from the two flags (see the mapping below).

| Input | Source |
| --- | --- |
| Phenols, 7 molecules (SDF) | ACC, `src/backend/tests/data/integrations/phenols.sdf` |
| 2HU4 with oseltamivir, hydrogens added (mmCIF; also converted to PDB with `gemmi convert`, which writes no `CONECT`) | ACC example `oseltamivir/2hu4.fw2.cif` |
| 8YAX, 40,738 atoms (mmCIF) | ACC example `pore/8yax.fw2.cif` |
| 1EK9 (mmCIF and converted PDB), 10AW | `tests/fixtures/corpus/cif/` |
| Representative small molecules, heme, fixed ions | `tests/fixtures/corpus/sdf/audit/representative.sdf`, `tests/fixtures/corpus/sdf/heme/ideal.sdf`, `tests/fixtures/synthetic/cif/fixed_ions.cif` |
| 1CBN (alternate locations, ligand), 1D3Z (10 NMR models), as mmCIF and PDB | downloaded from `files.rcsb.org` |

For rendering, the real ACC results page ran locally (`vite` dev server) against a mock backend that
serves only the method catalog and the per-molecule mmCIF files, in headless Chrome. For each file the
test recorded what Mol* built (atoms, bonds and orders, entity types, components and representations,
charge rows, atoms without a charge), ACC's own consistency check, and screenshots of the three view
types.

Not tested: malformed input records, calculations that fail numerically, MOL2 input, and modified
residues such as MSE in a PDB file downloaded from the PDB (the 1EK9 PDB file was converted by Gemmi).

## API mapping

| ChargeFW2 | ChargeFW 3 |
| --- | --- |
| `Molecules(file, ...)` | `chargefw.io.read(path, format=..., selection=..., conformers="first")`; `format` from the file extension; `selection` and `conformers` only for PDB and mmCIF |
| `read_hetatm=True, ignore_water=False` | `selection="all"` |
| `read_hetatm=True, ignore_water=True` | `selection="polymers-and-ligands"` |
| `read_hetatm=False` | `selection="polymers"` (not the same atoms, see below) |
| `permissive_types` | `parameter_matching="permissive"` or `"strict"` |
| `molecules.info()` | counts from `MoleculeCollection`: `len()`, `atom_count`, `atomic_numbers` |
| `get_available_methods()` | `chargefw.methods` |
| `get_available_parameters(method)` | `chargefw.methods[method].parameter_sets` |
| `get_suitable_methods(molecules)` | `chargefw.assess(molecules, method_scope="all", parameter_matching=...).plans`; keep one plan per (method, parameter set), because the list also has the cutoff and cover variants of each |
| `get_best_parameters(molecules, method, ...)` | `chargefw.assess(molecules, method=method, ...).default_plan.parameter_set` |
| `calculate_charges(...)` | `chargefw.calculate(molecules, method=..., parameter_set=..., parameter_matching=..., threads=1)`; charges are in `result.assignments`, in input order |
| written `.pqr`, `.mol2` | `chargefw.io.write(path, result, format="pqr" / "mol2")` |
| written `.txt` | no writer; ACC writes it from the assignments |
| written `.fw2.cif` | see [Viewer files](#viewer-files) |

Metadata ACC's models need:

- `Method`: `id` → `internal_name`, `name`, `full_name`, `publication`. `type` is `"3D"` when
  `requires_coordinates`, else `"2D"`; ChargeFW2 has `"other"` for `formal` and `dummy`.
  `has_parameters` is `len(method.parameter_sets) > 0`.
- `Parameters`: `id` → `internal_name`, `name` → `full_name`, `publication`, `notes`, `method.id`.

`chargefw.calculate` releases the GIL: four calculations in four threads took 2.2 s where one took
1.9 s. A full calculation uses one thread whatever `threads` says, so `threads=1` only keeps cutoff
and cover calculations from competing with ACC's own pool.

## Differences

### Catalogs

- The 22 method IDs are the same. Full names differ in capitalization and spelling for the SQE
  family; PEOE and MPEOE are "… of Orbital Electronegativity" (ChargeFW2 says "Atomic"); `formal` is
  "Formal" / "Formal atomic charges" and `dummy` is "Dummy method" / "Dummy zero charges".
- Parameter-set IDs lost the two-digit priority: `SQEqp_10_Schindler2021_CCD_gen` is
  `SQEqp_Schindler2021_CCD_gen`. `KCM_00_original` is `KCM_original_initial`, with the DOI corrected
  to `10.1002/jcc.20892`, and `KCM_original_fitted` is new. Names, publications, and notes of the
  other 51 sets are identical apart from American spelling in the notes.
- Method order is the same except that TSEF now precedes DENR.
- ChargeFW2 lists parameter sets derived from proteins first when the input contains a protein, so
  EEM on a protein starts with `EEM_Ionescu2013_npa_gas`. ChargeFW has one order for every input and
  starts with `EEM_NEEMP_ccd2016_npa`. The Python interface does not say which sets are protein sets.

Where ACC shows a parameter-set ID, checked in the frontend code:

| Place | Shown |
| --- | --- |
| Parameter selector on the setup page | the name, such as "Schindler 2021 (CCD_gen)" |
| "Parameters:" line on the results page | the name; the raw ID when the API does not know the ID |
| Charge Set selector on the results page | the raw string from the mmCIF file, `sqeqp/SQEqp_10_Schindler2021_CCD_gen` |
| Calculation history | the raw IDs, `sqeqp (SQEqp_10_Schindler2021_CCD_gen)` |
| Public API and its documentation in `src/frontend/acc/src/assets/wiki.md` | the raw IDs |

Download file names were not checked.

### Charges of small molecules

On the phenols and the representative molecules, every method agrees with ChargeFW2 within 1e-4 e
except three that ChargeFW corrected:

| Method | Largest difference |
| --- | --- |
| QEq | 0.36 e |
| Charge2 | 0.09 e |
| MGC | 0.09 e (heme) |

### Charges of structures

EQeq, which uses no bonds, agrees exactly. Methods that use bonds differ because ChargeFW2 reads
structures wrongly in two ways:

- Residues that share a number and differ in insertion code are merged, and their bonds are lost. In
  2HU4 (169/169A, 345/345A, 412/412A–D) the carbonyl carbons of these residues get −0.08 e from
  ChargeFW2 and about +0.65 e from ChargeFW with SQE+qp.
- Ligand hydrogens whose names differ from the CCD stay unbonded and get 0.040 e each; in the
  prepared 2HU4 that is most hydrogens of the ligand. ChargeFW takes the ligand's bonds from the
  file's `_chem_comp_bond`.

For 2HU4 the largest differences are 0.87 e (SQE+qp), 0.03–0.21 e (the EEM sets), and 0.53 e (QEq,
which also includes the method correction). For 10AW, which has no hydrogens, SQE+qp differs by at
most 0.09 e.

### Applicability

- ChargeFW rejects SQE, MGC, DENR, VEEM, PEOE, MPEOE, GDAC, and KCM for a molecule in which a bonded
  part has a nonzero total formal charge, and TSEF when the atoms are not one bonded part. ChargeFW2
  calculated them, with a wrong total charge. 2HU4 goes from 36 (method, parameter set) pairs to 26,
  the representative molecules from 21 to 17.
- Above 20,000 atoms ChargeFW also offers SQE+qp and SQEq0 with a cutoff, and SQE+qp is the default;
  ChargeFW2's default for 8YAX was EEM. Both use a 12 Å cutoff automatically.
- Heme: ChargeFW accepts `EEM_Svob2007_cmet2` and `EEM_Svob2007_hm2`, with strict matching too, and
  EEM becomes the default instead of EQeq. Both sets have parameters for every atom type of heme,
  iron with a bond included, so ChargeFW is right; why ChargeFW2 reports "parameters don't cover the
  molecule" was not found.
- A structure with atoms at identical coordinates gets no 3D method; see the first ChargeFW issue.

### Reading structures

- `selection="polymers"` goes by the entities of the file, not by the `HETATM` flag, so modified
  residues stay. 1EK9 keeps its 15 MSE residues: 9,918 atoms instead of 9,798, and only 6 pairs remain
  applicable instead of 52, because few sets have selenium.
- `polymers-and-ligands` leaves out `HOH`, `WAT`, `H2O`, and `DOD`, not only `HOH`.
- Alternate locations: ChargeFW reads the first one in the file for each atom name of a residue;
  ChargeFW2 reads atoms without one or with `A`. 1CBN gives 673 atoms instead of 640; see the first
  ChargeFW issue.
- ChargeFW reads all models as conformers by default; `conformers="first"` gives ChargeFW2's
  behavior.
- Bonds of ligands and modified residues come from the file. A PDB file without `CONECT` leaves them
  unbonded with an `unbonded_residue_atoms` warning, and the calculation still runs: in 2HU4 converted
  to PDB, each ligand atom keeps its formal charge, the amine nitrogen +1.000 e. The PDB file of 1CBN
  from the PDB does have `CONECT` records for its ligand, so the gap concerns PDB files written by
  tools that drop them.
- The same warning names stray hydrogens that preparation tools leave in the input. In 8YAX, `H1` of
  LEU A 1195 is 1.08 Å from the carbonyl carbon of the previous residue and 2.69 Å from its own
  nitrogen. ChargeFW is right to report them. It makes no peptide bond across a chain gap; this was
  checked on 2HU4, 8YAX, and 1EK9.
- A PDB molecule is named by its file name, an mmCIF molecule by its data block. ACC stores uploads
  as `<sha256>_<name>`, so the hash would become the name of a PDB molecule.
- Names are kept as written: `2,4-dinitrophenol`, where ChargeFW2 has `2_4_dinitrophenol`.

### Output files

- TXT: no writer.
- MOL2: SYBYL atom types (`O.3`) instead of element symbols, atom names `O1`, charges with four
  decimals, formal charges in `UNITY_ATOM_ATTR`, coordinates in the shortest exact form, which gives
  `8e-04` for 0.0008.
- PQR: two `REMARK` lines, charges and radii with four decimals.
- mmCIF: one file with one data block per molecule, charges dictionary 1.1 (`method` and
  `parameter_set` are separate items, plus `software_name` and `software_version`), bonds with
  orders in `_chem_comp_bond` and `_struct_conn`. `chargefw.io.write` builds the block from the
  molecule, so nothing else of an uploaded structure file is in it.
  `chargefw.io.gemmi.attach_charges(document, result)` adds the charge categories to the unchanged
  Gemmi document instead; it accepts only results calculated from that document
  (`chargefw.io.gemmi.from_document`).

### Speed

Fastest of three runs on a Ryzen 7 9700X, in seconds. ChargeFW2 in the ACC image is built with
`-march=native` and OpenMP, which lets Eigen parallelize its dense products. ChargeFW sets
`EIGEN_DONT_PARALLELIZE`, and its wheels target x86-64-v3.

| Calculation | ChargeFW2, 16 threads | ChargeFW2, 1 thread | rc5 wheel | rc5, local build, x86-64-v3 | rc5, local build, native |
| --- | --- | --- | --- | --- | --- |
| EEM, 2HU4, 5,805 atoms | 0.79 | 1.37 | 1.88 | 2.02 | 1.05 |
| EQeq, 1EK9, 11,426 atoms | 4.14 | 10.67 | 14.54 | 15.96 | 8.15 |
| SQE+qp, 2HU4 | 0.90 | 1.67 | 1.20 | 1.21 | 0.77 |
| QEq with cutoff, 8YAX, 40,738 atoms | 25.9 | – | 8.3 | – | – |

The implementation is not slower: built for the host CPU, ChargeFW beats single-threaded ChargeFW2
in all three full calculations. The wheel loses on EEM and EQeq because of its instruction set (AVX2
instead of this CPU's AVX-512). ChargeFW2's remaining lead is thread parallelism, which helps one
calculation on an idle machine and less on a server that runs several at once.

## Viewer files

Tested variants of the per-molecule file ACC serves to Mol*, against the file ChargeFW2 produces:

| Variant | Result |
| --- | --- |
| Structure written by ChargeFW (`io.write`, `io.gemmi.to_document`) | Worse. 8YAX loses its biological assembly: 40,738 atoms shown instead of the 244,428 of the pore. From PDB input, water and ligands become part of the polymer component. No secondary structure from the file |
| Charges attached to the unchanged upload | Fails when some atoms were not calculated (other alternate locations, water left out): ACC's check in `viewer.ts` throws "Atom count (777) does not match charge count (673)" |
| Attached, then the first model's atom rows without a charge removed | Same components, representations, and bonds as ChargeFW2 for 2HU4, 8YAX, 1EK9 with and without water, 1D3Z, and 1CBN |
| PDB converted to mmCIF by Gemmi, `_chem_comp` kept | Worse. Mol* shows everything as spacefill |
| PDB converted to mmCIF by Gemmi, `_chem_comp` removed | Same as ChargeFW2 |
| Small molecules written by ChargeFW | Same picture, and Mol* gets the bond orders (ChargeFW2's `_chem_comp_bond` rows use numeric atom IDs that match no atom name, so Mol* guesses the bonds) |
| Charge metadata as ChargeFW writes it (dictionary 1.1) | Worse. The Charge Set selector shows only `sqeqp`, and the "Parameters:" line is empty, because ACC splits `method` at `/` and the Mol* extension reads only `method` |

Details that matter:

- The Mol* extension matches `_sb_ncbr_partial_atomic_charges.atom_id` to `_atom_site.id` as
  integers. When the upload is kept, the charge rows must carry the file's atom IDs, not position plus
  one. `attach_charges` does that.
- With `_chem_comp` kept in an mmCIF from the PDB, Mol* additionally shows modified residues as
  ball-and-stick (120 MSE atoms in 1EK9). ChargeFW2 removes the category, so this is a choice.
- For a PDB upload, calculating from the Gemmi-converted document loses the `CONECT` bonds (1CBN:
  678 bonds instead of 680), and `attach_charges` refuses a result read from the PDB file. Until
  ChargeFW handles this (issue 3 below), ACC has to calculate from the PDB file and write the charge
  rows into the converted document itself. Both list the atoms in the same order.
- Where charges changed, the colors change with them: the scale maximum for 2HU4 goes from 1.10 to
  1.36.

## ChargeFW work

In the order proposed to the maintainer; none is started.

1. **Alternative residues at one position (bug).** 1CBN has PRO or SER at position 22 and LEU or ILE
   at 25, as alternate locations A and C. ChargeFW picks an alternate location for each atom name of
   a residue, so it reads both residues of each pair. Their backbone atoms coincide, and every method
   that uses coordinates is rejected with `invalid_geometry`, "identical coordinates"; the default
   becomes PEOE. ChargeFW2 keeps alternative A and calculates SQE+qp. The choice has to be made for
   the residue position.
2. **mmCIF output that keeps an mmCIF input.** Writing mmCIF for a result read from mmCIF should give
   the source data block with the charges, as `attach_charges` does, without the atom rows that were
   not calculated. Today the command line's `--output-mmcif` gives the rebuilt block, with the losses
   in the table above, for every user who opens the result in Mol*. `attach_charges` needs the same
   treatment of uncalculated atoms. This changes what the writer produces and needs a design sketch
   first.
3. **mmCIF output from PDB input.** The written block should have the entities, secondary structure,
   and assembly that Gemmi derives from the PDB file, as ChargeFW2's has, without `_chem_comp`, and
   with the bonds ChargeFW used. It belongs to the same design as item 2.
4. **MOL2 coordinates** in fixed notation, without exponents.

Decided against: CCD templates for all ligands (large data for a gap that the files from the PDB do
not have); sanitized names, TXT output, and per-molecule files (ACC's conventions); calculating
bond-based methods on charged molecules; leaving out molecules that fail; selection by `HETATM` flag.

Open questions for ChargeFW:

- Whether full calculations should use several threads. The `threads` option exists and does not
  affect them.
- Whether ACC should build ChargeFW from the sdist with
  `-C cmake.define.CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON` in its image, as it builds ChargeFW2 with
  `-march=native` today. On the test machine that makes full calculations 1.6–2 times faster.

Outside ChargeFW: the Mol* partial-charges extension should read `parameter_set` (dictionary 1.1).
Until it does, ACC must keep writing `method/parameters` into `method`.

## ACC work

Needed whatever is decided:

- Replace `integrations/chargefw2/` and the `chargefw2` imports in `services/chargefw2.py` with the
  calls in [API mapping](#api-mapping). `Molecules` objects are passed around between the services;
  a `MoleculeCollection` takes their place.
- Read structures with `conformers="first"`.
- Keep ChargeFW2's molecule names (sanitized, at most 80 characters, unique, structures named by
  entry ID) for dictionary keys, file names, and the `molecule` parameter.
- Write the TXT file: the name in upper case on one line, the charges with five decimals on the next.
- Produce the viewer files by the variants that match ChargeFW2 in [Viewer files](#viewer-files), one
  file per molecule, and keep writing the charge metadata as today.
- Stop reusing calculations stored before the switch.
- Dependencies: remove the `ChargeFW2` submodule and the `chargefw2-base:acc` image from
  `src/backend/Dockerfile`, `.github/workflows/deploy-logic.yml`, and the deployment documentation;
  add `chargefw` to `pyproject.toml` and `uv.lock`. ACC already depends on Gemmi 0.7.
- Regenerate `src/backend/app/examples/`; update `tests/services/test_chargefw2_service.py` and
  `tests/data/integrations/chargefw2_data.py`.
- Update the parameter-set IDs in `wiki.md` and the ChargeFW2 links in `wiki.md` and
  `components/home/sections/about.tsx`.

## Decisions

Decided by the maintainer:

- The priority numbers in the parameter-set IDs were a mistake; the new IDs are not a problem.
- Mol* rendering must not get worse.
- Corrections of ChargeFW2 errors are acceptable.

Open, with the recommendation given:

| Question | Recommendation |
| --- | --- |
| Should ACC's API still accept the old parameter-set IDs, and what happens to the IDs stored in the database? | Accept the old IDs as aliases; rename the stored ones |
| PDB files without `CONECT`: ligands unbonded | Show ChargeFW's `unbonded_residue_atoms` warning in ACC, as it shows the missing-hydrogens warning |
| Parameter-set order for proteins (EEM starts with a general set, not a protein set) | Accept ChargeFW's order |
| SQE+qp as the default above 20,000 atoms | Accept |
| `KCM_original_fitted` in the lists | Accept, or hide in ACC |
| "Read HETATM" off keeps modified residues | Accept |
| Fewer methods for charged molecules | Accept |
| A calculation fails when one molecule fails numerically | Accept |
| MOL2 and PQR downloads in ChargeFW's form | Accept |
| Users' stored calculation history from ChargeFW2 | Keep it viewable; do not reuse it for new requests |
| Keep `_chem_comp` in mmCIF uploads (modified residues shown as sticks) | Remove it, as ChargeFW2 does |
| ChargeFW issues 1–3 before 3.0.0, in another release candidate | Yes: 2 and 3 change the mmCIF output, which is easier before the stable release |
