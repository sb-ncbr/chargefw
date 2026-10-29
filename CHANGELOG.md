# Changelog

## Unreleased

### Command line

- A MOL2, mmCIF, PQR, or result JSON file that cannot be written no longer replaces an existing file of
  the same name. Before, a refused export, such as MOL2 for a molecule without coordinates, left the
  earlier file empty.
- A method option given as NaN or infinity, such as `--method-option denr.step=nan`, is refused as
  invalid input, exit status 2; it used to pass validation and fail later during the calculation.

### Python

- `chargefw.io.rdkit.from_mol()` refuses an RDKit molecule with a 2D conformer. Embed 3D coordinates,
  remove the conformers, or mark a real flat geometry with `Set3D(True)`.
- `molecule.import_diagnostics` lists the warnings a reader attached to an imported molecule, such as
  `partial_charges_ignored`, as `chargefw.ImportDiagnostic` values with a code, a message, and a line.
  Before, they were only in the result JSON.
- `chargefw.calculate()` with an assessed plan refuses `parameter_matching`, `execution`,
  `cutoff_threshold`, and `cover_threshold` also when they are passed with their default value, such as
  `parameter_matching="strict"`. Before, such arguments were silently ignored and the plan's own
  settings were used.

### C++

- Assessment refuses a custom method whose requirements list common parameters but no atom or bond
  parameters, with `std::invalid_argument`; declare such constants as method options. Common
  parameters belong in a parameter set only when they are fitted together with atom or bond
  parameters. Such a method used to pass assessment and then fail during execution.
  `MethodRequirements::requires_parameter_classification()` states which methods need atom or bond
  parameters.

### Output

- Writing mmCIF for many molecules with the same name, such as an SDF library whose records share a
  title, no longer slows down with the cube of their number: 4,000 such molecules took about 100
  seconds and now take half a second. Result JSON and MOL2 output of large collections is also
  faster.
- ChargeFW writes PQR files for electrostatics with APBS: `--output-pqr` on the command line,
  `format="pqr"` in Python, and `PqrWriter` in C++. PQR is available for one structure read from PDB or
  mmCIF, with one conformer; each atom gets its element's van der Waals radius. ChargeFW refuses
  structures that APBS would misread, such as chain IDs other than a single letter.
  [File formats](https://chargefw.readthedocs.io/en/latest/reference/formats.html#pqr) lists the rules.
- mmCIF output contains the bonds used in the calculation, with their orders, so viewers that read
  mmCIF bond data can show double and triple bonds. Bonds within a residue are in `_chem_comp_bond`;
  bonds between residues, such as peptide links and disulfides, are in `_struct_conn`. ChargeFW refuses
  to write a structure whose bonds the residue definitions cannot express, for example two residues of
  the same name with different bonds.
  [File formats](https://chargefw.readthedocs.io/en/latest/reference/formats.html#mmcif) lists the rules.
- An mmCIF file can hold several charge sets for the same molecules: pass a list of results to
  `chargefw.io.dumps()`, `chargefw.io.write()`, `chargefw.io.gemmi.to_document()`, or
  `chargefw.io.gemmi.attach_charges()`. Two sets with the same method and parameter set are refused.
- A charge set calculated per conformer now has one ID in mmCIF output, covering all conformers. Before,
  each conformer got its own ID, so viewers listed the same method once per conformer.
- Attaching charges to an mmCIF document replaces its other charges dictionary declarations in
  `_audit_conform`, such as `mmcif_charges_v10.dic` from ACC or ChargeFW2 files, and repeated ones;
  declarations of other dictionaries stay.

### Input

- An SDF or MOL record that declares a 2D drawing, such as a PubChem 2D download, is read without
  coordinates and gets a `coordinates_2d_ignored` warning. Methods that need 3D geometry no longer
  calculate charges from the flat drawing; the others still apply.
- MOL and SDF atoms no longer get their element symbol as a name, since molfiles have no atom names;
  `atom_names` in Python are empty for them. MOL2 and mmCIF output name these atoms from the element
  and atom position (`O1`, `O2`, `N4`), where they used to show bare symbols (`O`, `O`) or numbered
  repeats (`O`, `O2`).
- `chargefw.io.gemmi.from_structure()` uses template bonds by default and refuses `bonds="explicit"`
  and `"hybrid"`. Converting a `gemmi.Structure` loses the file's ligand bond orders, connection bond
  orders, and PDB `CONECT` records, so those bonds were lost or became single without a warning. With
  the new default, links that only the file records, such as disulfides and covalent ligands, are no
  longer bonded. For bonds from the file, use `from_document()` for mmCIF and `chargefw.io.read()` for
  PDB.
- PDB and mmCIF connections to an atom of a symmetry-related copy (`LINK`, `SSBOND`, or
  `_struct_conn` partners with different symmetry codes) are left out with a
  `symmetry_connections_ignored` warning. Before, they became bonds to the atom of the same name in the
  file, which is a different atom.
- An mmCIF block that contains partial charges, for example from an earlier calculation, gives its
  record a `partial_charges_ignored` warning, as MOL2 input does; the charges are not used.
- mmCIF bonds between residues keep their order from `_struct_conn.pdbx_value_order`; they were all
  read as single bonds before. A connection without an order is still a single bond.

## 3.0.0a4

### Installation

- ChargeFW requires Python 3.12 or newer. One wheel serves all of them, including Python versions
  released later. On Python 3.11, pip installs an older ChargeFW release instead;
  [Installation](https://chargefw.readthedocs.io/en/latest/getting-started/installation.html) shows how
  to create an environment with a newer Python.
- `pip install chargefw` also installs the `chargefw` command line.
- Wheels are also published for Linux aarch64, macOS 14 or newer on Apple silicon, and Windows
  x86-64. On Windows, the command line accepts file names that the system code page cannot represent.

### Distribution

- The wheel contains the ChargeFW library and Gemmi in its Python extension and command line, and no
  longer ships separate `libchargefw_core.so` and `libgemmi_cpp.so` libraries. It is about a quarter
  smaller.

### Building from source

- The new CMake option `CHARGEFW_BUILD_LIBRARY`, on by default, builds and installs the shared C++
  library. The library and its installation are supported on Linux; on Windows, CMake stops with an
  explanation and points to the Python package instead.
- ChargeFW compiles with Clang and libc++ 19 or newer, including Apple Clang from Xcode 16.4 and 26.1.

### Python

- The type aliases in `chargefw.core` and `chargefw.io`, such as `chargefw.io.InputFormat`, are defined
  with the `type` statement. Type checkers treat them as before; at run time, `typing.get_args()` needs
  the alias's `__value__`, for example `typing.get_args(chargefw.io.InputFormat.__value__)`.

## 3.0.0a3

### Documentation

- The [documentation](https://chargefw.readthedocs.io/en/latest/) is reorganized and rewritten for
  users of the Python package, the command line, and the C++ library. Pages have new addresses, so
  links to the old ones, such as `PYTHON.html`, `FORMATS.html`, or `recipes/`, no longer work.
  [Differences from ChargeFW2](https://chargefw.readthedocs.io/en/latest/about/chargefw2.html)
  describes what changed from ChargeFW2 and how to move over.

### Output

- MOL2 and mmCIF files write charges rounded to four decimal places, and coordinates with the digits
  they were read with instead of 17 significant digits. Result JSON keeps full precision.

### Listings

- `chargefw methods`, `chargefw.methods`, and the C++ method registry list methods by descending
  priority, then ID. `chargefw parameters`, `chargefw.parameter_sets`, and the C++ parameter loaders list
  parameter sets grouped by method, each method's sets by descending priority, then ID.

## 3.0.0a2

### Building from source

- Native builds use installed CLI11, nlohmann/json, Eigen, nanoflann, oneTBB, and Gemmi packages by
  default and report any that are missing at configure time. The supported versions are listed in
  the [native library documentation](https://chargefw.readthedocs.io/en/latest/interfaces/cpp.html#build-and-link);
  the minimums match Debian 13. Set
  `CHARGEFW_FETCH_DEPENDENCIES=ON` to download pinned sources instead, as the developer presets and
  the Python package build do. This option replaces `CHARGEFW_USE_SYSTEM_DEPENDENCIES`; builds that
  set it to `ON` should use `CHARGEFW_FETCH_DEPENDENCIES=OFF`.
- Downloaded dependency sources are verified against SHA-256 checksums.
- The installed CMake package requires the exact nlohmann/json and Gemmi versions that ChargeFW was
  built against.
- Installed binaries find dependency libraries installed outside the system library paths without
  `LD_LIBRARY_PATH`.

### Distribution

- Wheels, native installations, and the CLI image include the license notices of the third-party
  code in their binaries.

### Fixes

- mmCIF input accepts tokens longer than 4096 bytes; entries with long values, such as 3J3Q, previously
  failed to parse.
- Cover execution of SQE+qp conserves the total charge on very large structures; on 3J3Q (2.44 million
  atoms), rounding had left the total at 1.1e-4 instead of 0.

## 3.0.0a1

First pre-release of the ChargeFW 3 series, the successor of
[ChargeFW2](https://github.com/sb-ncbr/ChargeFW2). The version number continues ChargeFW2's lineage;
the `chargefw` Python package moves from 0.1.7 directly to 3.0.0a1. Install it with
`pip install --pre chargefw`.

ChargeFW 3 is a new implementation. Its highlights are:

- one native engine behind a C++23 library, a command-line application, and a Python package for
  CPython 3.11 to 3.14 on Linux x86-64;
- 22 methods and the bundled parameter sets, with deterministic method and parameter-set selection and
  structured reasons for rejected candidates;
- explicit molecular input: topology, protonation, formal charges, and coordinates are never changed
  silently, and missing prerequisites are reported as failures;
- full calculations and explicit cutoff and cover approximations, recorded in result provenance;
- SDF, MOL2, PDB, mmCIF, and molecule JSON input; result JSON, MOL2, and mmCIF output; source atom,
  molecule, and conformer order preserved in every result;
- optional Gemmi and RDKit integration in Python.

Documentation is published at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/).

### Changes since 0.1.7

- Parameter-set keys are exposed in Python (`atom_keys` and `bond_keys`) and in the CLI
  (`chargefw parameters ID --keys`).
- The bundled periodic table is exposed in Python as `chargefw.periodic_table`.
- New recipes calculate dipole moments and attach charges to RDKit SDF records; the file calculation
  recipe can export CSV, and the fixed-ion and molecule-inspection recipes report more detail.
- The repository moved to [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw).
