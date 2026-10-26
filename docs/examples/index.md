# Examples

Each example is a complete Python script that you can download and run. They use the
[Python package](../interfaces/python.md) and show its typical uses:

| Example | What you get |
| --- | --- |
| [Calculate a molecular file](calculate_file.md) | Charges of every atom of a file, printed or as CSV |
| [Inspect imported molecules](inspect_molecules.md) | What ChargeFW reads from a file, and how parameter sets see its atoms |
| [Charge a Gemmi document](charge_gemmi_document.md) | An mmCIF file with charges added; for mmCIF input, all its other data is kept |
| [Calculate with a fixed ion](calculate_with_fixed_ions.md) | Charges of a structure with ions, and the charges around each ion |
| [Analyze conformer-dependent charges](analyze_rdkit_conformer_charges.md) | Which atoms' charges change most between RDKit conformers |
| [Calculate dipole moments](calculate_dipole_moments.md) | Dipole moments from the calculated charges |
| [Write charges to RDKit molecules](charge_rdkit_sdf.md) | An SDF file with charges stored as RDKit atom properties |
| [Compare bundled parameter sets](compare_parameter_sets.md) | How much the charges differ between the parameter sets of one method |

```{toctree}
:hidden:

calculate_file
inspect_molecules
charge_gemmi_document
calculate_with_fixed_ions
analyze_rdkit_conformer_charges
calculate_dipole_moments
charge_rdkit_sdf
compare_parameter_sets
```
