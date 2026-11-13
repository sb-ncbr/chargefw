#include "bindings.h"
#include "native_execution_result.h"
#include "native_input_metadata.h"

#include "io/gemmi/mmcif_attachment.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/io/charge_result.h>
#include <chargefw/io/conformer_selection.h>
#include <chargefw/io/json_input.h>
#include <chargefw/io/json_output.h>
#include <chargefw/io/mmcif_input.h>
#include <chargefw/io/mmcif_output.h>
#include <chargefw/io/mol2_input.h>
#include <chargefw/io/mol2_output.h>
#include <chargefw/io/mol_input.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/io/pdb_input.h>
#include <chargefw/io/pqr_output.h>
#include <chargefw/io/sdf_input.h>
#include <chargefw/io/structure_options.h>

#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <gemmi/cif.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

struct MoleculePayload {
    std::vector<int> atomic_numbers;
    std::vector<int> formal_charges;
    std::vector<std::array<std::size_t, 3>> bonds;
    std::vector<std::vector<std::array<double, 3>>> coordinates;
    std::string name;
    std::vector<std::string> atom_names;
    std::vector<std::string> conformer_names;
    io::MoleculeRecordIdentity identity;
    std::optional<io::MoleculeImportMetadata> import_metadata;
    std::vector<io::ImportIssue> import_warnings;
    std::shared_ptr<NativeInputMetadata> native_input_metadata;
};

[[nodiscard]] auto atom_references(const std::vector<io::SourceAtomReference>& references)
    -> nb::list {
    auto result = nb::list{};
    for (const auto& reference : references) {
        if (reference.structural_labels == nullptr) {
            result.append(nb::make_tuple(reference.position, reference.id, nb::none()));
            continue;
        }
        const auto hierarchy = [](const io::SourceHierarchyLabels& labels) {
            auto value = nb::dict{};
            value["atom"] = labels.atom;
            value["residue"] = labels.residue;
            value["chain"] = labels.chain;
            value["sequence"] = labels.sequence;
            return value;
        };
        const auto& labels = *reference.structural_labels;
        auto encoded = nb::dict{};
        encoded["author"] = hierarchy(labels.author);
        encoded["label"] = hierarchy(labels.label);
        encoded["entity"] = labels.entity;
        encoded["insertion_code"] = labels.insertion_code;
        encoded["alternate_location"] = labels.alternate_location;
        encoded["segment"] = labels.segment;
        result.append(nb::make_tuple(reference.position, reference.id, std::move(encoded)));
    }
    return result;
}

[[nodiscard]] auto import_metadata(const io::MoleculeImportMetadata& metadata) -> nb::dict {
    auto conformers = nb::list{};
    for (const auto& conformer : metadata.conformers) {
        conformers.append(
            nb::make_tuple(conformer.position, conformer.id, atom_references(conformer.sites)));
    }
    auto result = nb::dict{};
    auto components = nb::list{};
    for (const auto& component : metadata.components) {
        components.append(nb::make_tuple(component.component_id, component.atom_indices));
    }
    result["format"] = std::string{io::to_string(metadata.format)};
    result["atoms"] = atom_references(metadata.atoms);
    result["conformers"] = std::move(conformers);
    auto policy = nb::dict{};
    policy["selection"] = metadata.policy.selection;
    policy["alternate_locations"] = metadata.policy.alternate_locations;
    policy["conformers"] = metadata.policy.conformers;
    policy["bonds"] = metadata.policy.bonds;
    result["policy"] = std::move(policy);
    result["source_connectivity"] = std::string{io::to_string(metadata.source_connectivity)};
    result["components"] = std::move(components);
    return result;
}

[[nodiscard]] auto portable_id(const io::PortableId& id) -> nb::object {
    if (id.empty()) {
        return nb::none();
    }
    return std::visit([](const auto& value) { return nb::cast(value); }, *id.value());
}

auto make_payload(io::ImportedMoleculeRecord record) -> MoleculePayload {
    MoleculePayload result;
    const auto& molecule = record.molecule;
    result.atomic_numbers.reserve(molecule.atom_count());
    result.formal_charges.reserve(molecule.atom_count());
    result.atom_names.reserve(molecule.atom_count());
    for (const auto& atom : molecule.atoms()) {
        result.atomic_numbers.push_back(atom.atomic_number());
        result.formal_charges.push_back(atom.formal_charge());
        result.atom_names.emplace_back(atom.name());
    }

    result.bonds.reserve(molecule.bond_count());
    for (const auto& bond : molecule.bonds()) {
        result.bonds.push_back({bond.first_atom_index(), bond.second_atom_index(),
                                core::bond_order_value(bond.order())});
    }

    result.coordinates.reserve(molecule.conformer_count());
    result.conformer_names.reserve(molecule.conformer_count());
    for (const auto& conformer : molecule.conformers()) {
        auto positions = std::vector<std::array<double, 3>>{};
        positions.reserve(conformer.size());
        for (const auto& position : conformer.positions()) {
            positions.push_back({position.x, position.y, position.z});
        }
        result.coordinates.push_back(std::move(positions));
        result.conformer_names.emplace_back(conformer.name());
    }

    result.name = molecule.name();
    result.import_warnings = record.import_warnings;
    result.native_input_metadata = std::make_shared<NativeInputMetadata>(
        record.identity, record.import_warnings, record.import_metadata);
    result.identity = std::move(record.identity);
    result.import_metadata = std::move(record.import_metadata);
    return result;
}

auto as_python(const MoleculePayload& payload) -> nb::dict {
    nb::dict result;
    result["atomic_numbers"] = nb::cast(payload.atomic_numbers);
    result["formal_charges"] = nb::cast(payload.formal_charges);
    result["bonds"] = nb::cast(payload.bonds);
    result["coordinates"] = nb::cast(payload.coordinates);
    result["name"] = payload.name;
    result["atom_names"] = nb::cast(payload.atom_names);
    result["conformer_names"] = nb::cast(payload.conformer_names);
    result["source_name"] = payload.identity.source_name;
    result["record_index"] = payload.identity.record_index;
    result["record_id"] = portable_id(payload.identity.record_id);
    nb::list import_warnings;
    for (const auto& issue : payload.import_warnings) {
        import_warnings.append(nb::make_tuple(
            std::string{io::to_string(issue.kind)}, issue.message,
            issue.line.has_value() ? nb::cast(*issue.line) : nb::none(),
            nb::cast(issue.location.atom_indices), nb::cast(issue.location.bond_indices)));
    }
    result["import_warnings"] = import_warnings;
    result["native_input_metadata"] = payload.native_input_metadata;
    if (payload.import_metadata.has_value()) {
        result["import_metadata"] = import_metadata(*payload.import_metadata);
    } else {
        result["import_metadata"] = nb::none();
    }
    return result;
}

auto as_python(const std::vector<MoleculePayload>& payloads) -> nb::list {
    nb::list result;
    for (const auto& payload : payloads) {
        result.append(as_python(payload));
    }
    return result;
}

template <typename Reader> auto read_all(Reader& reader) -> std::vector<MoleculePayload> {
    auto payloads = std::vector<MoleculePayload>{};
    while (auto record = reader.next()) {
        payloads.push_back(make_payload(std::move(*record)));
    }
    return payloads;
}

auto parse(std::string contents, std::string source_name, const std::string& format,
           const std::string& selection, const std::string& bonds, const std::string& conformers)
    -> nb::list {
    auto payloads = std::vector<MoleculePayload>{};
    {
        nb::gil_scoped_release release;
        std::istringstream input{std::move(contents)};
        const auto structural_options = [&] {
            return io::StructureOptions{.selection = io::record_selection_from_string(selection),
                                        .bonds = io::bond_strategy_from_string(bonds),
                                        .conformers =
                                            io::conformer_selection_from_string(conformers)};
        };
        switch (io::molecular_source_format_from_string(format)) {
        case io::MolecularSourceFormat::mol: {
            auto reader = io::MolReader{input, std::move(source_name)};
            payloads = read_all(reader);
            break;
        }
        case io::MolecularSourceFormat::sdf: {
            auto reader = io::SdfReader{input, std::move(source_name)};
            payloads = read_all(reader);
            break;
        }
        case io::MolecularSourceFormat::mol2: {
            auto reader = io::Mol2Reader{input, std::move(source_name)};
            payloads = read_all(reader);
            break;
        }
        case io::MolecularSourceFormat::molecule_json: {
            auto reader = io::JsonReader{input, std::move(source_name),
                                         io::conformer_selection_from_string(conformers)};
            payloads = read_all(reader);
            break;
        }
        case io::MolecularSourceFormat::pdb: {
            auto reader = io::PdbReader{input, std::move(source_name), structural_options()};
            payloads = read_all(reader);
            break;
        }
        case io::MolecularSourceFormat::mmcif: {
            auto reader = io::MmcifReader{input, std::move(source_name), structural_options()};
            payloads = read_all(reader);
            break;
        }
        }
    }
    return as_python(payloads);
}

[[nodiscard]] auto result_pointers(const std::vector<NativeExecutionResult*>& native_results)
    -> std::vector<const io::ChargeCalculationResult*> {
    auto results = std::vector<const io::ChargeCalculationResult*>{};
    results.reserve(native_results.size());
    for (const auto* native_result : native_results) {
        results.push_back(std::addressof(native_result->result()));
    }
    return results;
}

auto dumps(const std::vector<NativeExecutionResult*>& native_results, const std::string& format)
    -> std::string {
    auto output = std::ostringstream{};
    if (format == "mmcif") {
        const auto results = result_pointers(native_results);
        io::MmcifWriter{output}.write(results);
        return output.str();
    }
    if (native_results.size() != 1) {
        throw std::invalid_argument{"output format " + format + " holds one calculation result"};
    }
    const auto& result = native_results.front()->result();
    if (format == "result-json") {
        io::JsonWriter{output}.write(result);
    } else if (format == "mol2") {
        io::Mol2Writer{output}.write(result);
    } else if (format == "pqr") {
        io::PqrWriter{output}.write(result);
    } else {
        throw std::invalid_argument{"unsupported calculation output format: " + format};
    }
    return output.str();
}

auto attach_mmcif(std::string contents, const std::vector<NativeExecutionResult*>& native_results,
                  const bool overwrite) -> std::string {
    auto output = std::ostringstream{};
    const auto document =
        ::gemmi::cif::read_memory(contents.data(), contents.size(), "<Gemmi document>");
    const auto results = result_pointers(native_results);
    io::gemmi::mmcif_output::write_attached(output, results, document, overwrite);
    return output.str();
}

} // namespace

void bind_io(nb::module_& module) {
    [[maybe_unused]] const auto native_input_metadata =
        nb::class_<NativeInputMetadata>(module, "_NativeInputMetadata");
    module.def("_parse", &parse, nb::arg("contents"), nb::arg("source_name"), nb::arg("format"),
               nb::arg("selection"), nb::arg("bonds"), nb::arg("conformers"));
    module.def("_dumps", &dumps, nb::arg("results"), nb::arg("format"),
               nb::call_guard<nb::gil_scoped_release>());
    module.def("_attach_mmcif", &attach_mmcif, nb::arg("contents"), nb::arg("results"),
               nb::arg("overwrite"), nb::call_guard<nb::gil_scoped_release>());
}

} // namespace chargefw::python
