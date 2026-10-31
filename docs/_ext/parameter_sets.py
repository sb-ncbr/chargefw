"""Generate the table of bundled parameter sets from data/parameters before Sphinx reads the sources."""

import json
import re
from collections import defaultdict
from pathlib import Path

ATOM_TYPES = {
    "plain": "element",
    "hbo": "element, highest bond order",
    "bonded": "element, bonded neighbors",
}
BOND_SYMBOLS = {"1": "–", "2": "=", "3": "≡"}
METADATA = re.compile(r"MethodMetadata\s+metadata\{(.*?)\};", re.DOTALL)


def string_field(block, field):
    # Adjacent literals are how clang-format wraps long strings.
    match = re.search(rf"\.{field}\s*=\s*((?:\"[^\"]*\"\s*)+)", block)
    return "".join(re.findall(r"\"([^\"]*)\"", match.group(1))) if match else None


def method_names(source_root):
    names = {}
    for header in sorted((source_root / "src" / "methods" / "builtin").glob("*.h")):
        for block in METADATA.findall(header.read_text(encoding="utf-8")):
            method_id = string_field(block, "id")
            names[method_id] = (string_field(block, "name"), string_field(block, "full_name"))
    return names


def bond_column(document):
    bonds = (document.get("bond") or {}).get("data", [])
    if not bonds:
        return ""
    if any(atom["classifier"] != "plain" for entry in bonds for atom in entry["key"]["atoms"]):
        return f"{len(bonds)} entries by atom types and bond order"
    listed = []
    for entry in bonds:
        first, second = (atom["element"] for atom in entry["key"]["atoms"])
        if first == "*" and second == "*":
            listed.append("any bond")
        else:
            listed.append(first + BOND_SYMBOLS.get(entry["key"]["bond"]["type"], "–") + second)
    if (
        all(entry["key"]["bond"]["classifier"] == "plain" for entry in bonds)
        and "any bond" not in listed
    ):
        return ", ".join(listed) + "; any bond order"
    return ", ".join(listed)


def table_row(document, with_bonds):
    metadata = document["metadata"]
    atoms = (document.get("atom") or {}).get("data", [])
    elements = sorted({entry["key"]["element"] for entry in atoms})
    classifiers = sorted({entry["key"]["classifier"] for entry in atoms})
    name = metadata.get("name", metadata["id"])
    if metadata.get("publication"):
        name = f"[{name}](https://doi.org/{metadata['publication']})"
    types = "; ".join(ATOM_TYPES[classifier] for classifier in classifiers)
    bonds = f" {bond_column(document)} |" if with_bonds else ""
    return (
        f"| {name} | `{metadata['id']}` | {', '.join(elements)} | {types} |{bonds} "
        f"{metadata.get('notes', '')} |"
    )


def render(source_root):
    names = method_names(source_root)
    by_method = defaultdict(list)
    for path in sorted((source_root / "data" / "parameters").glob("*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        by_method[document["metadata"]["method"]].append(document)

    lines = []
    for method_id in sorted(by_method, key=lambda key: names.get(key, (key,))[0].lower()):
        name, full_name = names.get(method_id, (method_id, None))
        lines += [f"## {name}", ""]
        if full_name and full_name != name:
            lines += [f"{full_name}, method `{method_id}`.", ""]
        else:
            lines += [f"Method `{method_id}`.", ""]
        with_bonds = any(
            (document.get("bond") or {}).get("data") for document in by_method[method_id]
        )
        if with_bonds:
            lines += [
                "| Parameter set | ID | Elements | Atom types by | Bonds | Notes |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        else:
            lines += [
                "| Parameter set | ID | Elements | Atom types by | Notes |",
                "| --- | --- | --- | --- | --- |",
            ]
        # The order in which automatic selection tries the sets.
        documents = sorted(
            by_method[method_id],
            key=lambda document: (
                -document["metadata"].get("priority", 0),
                document["metadata"]["id"],
            ),
        )
        lines += [table_row(document, with_bonds) for document in documents]
        lines.append("")
    return "\n".join(lines)


def generate(app):
    docs = Path(app.confdir)
    output = docs / "_generated" / "parameter-sets.md"
    output.parent.mkdir(exist_ok=True)
    content = render(docs.parent)
    if not output.exists() or output.read_text(encoding="utf-8") != content:
        output.write_text(content, encoding="utf-8")


def setup(app):
    app.connect("builder-inited", generate)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
