# Preparing input

ChargeFW calculates empirical partial atomic charges from an explicit molecular graph, formal charges,
and, when required by a method, Cartesian coordinates. The same native engine is used by the C++ API,
the CLI, and the Python package.

**Explicit molecular input.** ChargeFW does not silently change topology, protonation, formal charge,
hydrogen count, coordinates, method options, parameter sets, or execution mode.

The language-independent [molecular format reference](../reference/formats.md) defines supported inputs and outputs,
exactly what each reader imports, and how each writer preserves or generates molecular data.

ChargeFW is not a general chemistry-preparation toolkit. It does not provide SMILES parsing, arbitrary
bond perception, sanitization, protonation, hydrogen addition or removal, conformer generation, or
geometry optimization. Structural template bonding covers common amino acids, nucleotides, water, and
basic polymer links; it is not a complete Chemical Component Dictionary provider.
