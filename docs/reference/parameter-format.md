# Parameter-set format

ChargeFW's parameter sets are JSON files installed with ChargeFW; the
[parameter-set reference](parameter-sets.md) lists them. This page describes the file format and how a
parameter set matches atoms and bonds, for people who write or check parameter sets. The command line
uses the installed sets. Python and C++ can also load your own files, as described in [Methods and
parameter sets](../interfaces/python.md#methods-and-parameter-sets) for Python and for
[C++](../interfaces/cpp.md#methods-and-parameter-sets).

A parameter set holds the values for one method. The method determines which common, atom, and bond
parameter names must be present. Before a calculation, ChargeFW checks the required names and matches
the atoms and bonds that need parameters to entries in the set, as described below. Atoms kept as
[fixed ions](../guide/fixed-ions.md) are not matched.

## Document structure

```json
{
  "metadata": {
    "id": "example",
    "method": "sqeqp",
    "name": "Example",
    "publication": "10.0000/example",
    "notes": "Optional provenance",
    "priority": 0
  },
  "common": {
    "names": ["scale"],
    "values": [1.0]
  },
  "atom": {
    "names": ["electronegativity", "hardness", "width", "q0"],
    "data": [
      {
        "key": {"element": "C", "classifier": "hbo", "type": "1"},
        "values": [5.0, 10.0, 1.0, 0.0]
      }
    ]
  },
  "bond": {
    "names": ["kappa"],
    "data": [
      {
        "key": {
          "atoms": [
            {"element": "C", "classifier": "hbo", "type": "1"},
            {"element": "O", "classifier": "hbo", "type": "1"}
          ],
          "bond": {"classifier": "bo", "type": "1"}
        },
        "values": [0.5]
      }
    ]
  }
}
```

`metadata.id` and `metadata.method` are required. The non-empty ID names the parameter set, not the
file name, and no two loaded parameter sets may share it. `name`, `publication`, `notes`, and
`priority` are optional. `priority` is an integer from `0` to `65535`; automatic selection tries
sets with a higher priority first, and orders sets with equal priority by their ID.

`common`, `atom`, and `bond` are optional. Each has a `names` array that gives the order of its
`values`: every `values` array has one finite number for each name, in that order. A name may appear
only once in its section.

Members not described on this page are ignored.

## Atom parameters

Each `atom.data` entry has an atom key and a values array:

```json
{
  "key": {"element": "C", "classifier": "bonded", "type": "CNO"},
  "values": [2.6846, 1.571, 0.5771, 0.8167]
}
```

The key fields are all required:

- `element` is a supported element symbol, or `"*"` to match every element.
- `classifier` selects how ChargeFW derives the atom's `type`.
- `type` is the expected string produced by that classifier.

## Bond parameters

Each `bond.data` entry uses two atom keys and one bond key:

```json
{
  "key": {
    "atoms": [
      {"element": "C", "classifier": "bonded", "type": "CNO"},
      {"element": "N", "classifier": "bonded", "type": "CCH"}
    ],
    "bond": {"classifier": "bo", "type": "1"}
  },
  "values": [0.5]
}
```

`atoms` has exactly two atom keys, and a bond matches the entry with its atoms in either order.
`bond` needs `classifier` and `type`, which describe the bond itself.

Some bond parameters belong to one of the two atoms. ABEEM applies `D` to the atom with the higher
atomic number, or with the lower atom index when both are the same element, and `C` to the other
one, whatever the order of the keys.

## Classifiers

Atom keys support these classifiers:

| Classifier | Derived type | Example |
| --- | --- | --- |
| `plain` | Always `"*"` | `{"element": "O", "classifier": "plain", "type": "*"}` |
| `hbo` | The highest order of the atom's bonds, as a string; `"0"` for an atom without bonds | An atom with a double bond has type `"2"` |
| `bonded` | The element symbols of the bonded atoms, one per bond, in alphabetical order; `""` for an atom without bonds | Carbon bonded to C, N, and O has type `"CNO"`; a methyl carbon bonded to another carbon `"CHHH"` |

Bond keys support these classifiers:

| Classifier | Derived type | Example |
| --- | --- | --- |
| `plain` | Always `"*"` | `{"classifier": "plain", "type": "*"}` |
| `bo` | Bond order as `"1"`, `"2"`, or `"3"` | A double bond has type `"2"` |

A `plain` key must have the type `"*"`; any other type is an error when the file is loaded. `bonded`
uses only the elements of the bonded atoms, not bond orders or formal charges.

## Entry selection

ChargeFW goes through the entries in the order of `data` and uses the first one that matches. Put a
more specific key before a broader one to give it precedence. Repeated or overlapping keys are
allowed; a later entry that also matches is not used.

## Strict and permissive matching

Strict matching is the default. With permissive matching (`--permissive-types` on the command line,
`parameter_matching="permissive"` in Python), an atom or bond that no entry matches tries again with
its `hbo` or `bo` type lowered by one: an atom with a double bond, type `"2"`, can then use an entry
for `"1"`, and a triple bond an entry for `"2"`. An atom without bonds, type `"0"`, can use an entry
for `"1"`. `plain` and `bonded` keys match the same way in both modes.
