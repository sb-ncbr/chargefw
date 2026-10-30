#include "fixed_ions.h"

#include "component.h"
#include "core/diagnostic_description.h"

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace chargefw::adapters::detail {
namespace {

struct Ion {
    std::string_view component;
    int atomic_number;
    double charge;
    bool common = false;
};

// CCD snapshot records: https://files.rcsb.org/pub/pdb/data/monomers/components.cif.gz
// Verified on 2026-10-05: released monatomic ions with nonzero atom/component charge agreement.
constexpr std::array ions{
    Ion{"NA", 11, 1.0, true}, Ion{"K", 19, 1.0, true},   Ion{"MG", 12, 2.0, true},
    Ion{"CA", 20, 2.0, true}, Ion{"CL", 17, -1.0, true}, Ion{"ZN", 30, 2.0, true},
    Ion{"FE", 26, 3.0, true}, Ion{"FE2", 26, 2.0, true}, Ion{"0BE", 4, 2.0},
    Ion{"3CO", 27, 3.0},      Ion{"3NI", 28, 3.0},       Ion{"4MO", 42, 4.0},
    Ion{"4PU", 94, 4.0},      Ion{"4TI", 22, 4.0},       Ion{"6MO", 42, 6.0},
    Ion{"AG", 47, 1.0},       Ion{"AL", 13, 3.0},        Ion{"AM", 95, 3.0},
    Ion{"AU", 79, 1.0},       Ion{"AU3", 79, 3.0},       Ion{"BA", 56, 2.0},
    Ion{"BR", 35, -1.0},      Ion{"BS3", 83, 3.0},       Ion{"CD", 48, 2.0},
    Ion{"CE", 58, 3.0},       Ion{"CF", 98, 3.0},        Ion{"CO", 27, 2.0},
    Ion{"CR", 24, 3.0},       Ion{"CS", 55, 1.0},        Ion{"CU", 29, 2.0},
    Ion{"CU1", 29, 1.0},      Ion{"CU3", 29, 3.0},       Ion{"D8U", 1, 1.0},
    Ion{"DY", 66, 3.0},       Ion{"ER3", 68, 3.0},       Ion{"EU", 63, 2.0},
    Ion{"EU3", 63, 3.0},      Ion{"F", 9, -1.0},         Ion{"GA", 31, 3.0},
    Ion{"GD3", 64, 3.0},      Ion{"HG", 80, 2.0},        Ion{"HO3", 67, 3.0},
    Ion{"IN", 49, 3.0},       Ion{"IOD", 53, -1.0},      Ion{"IR", 77, 4.0},
    Ion{"IR3", 77, 3.0},      Ion{"LA", 57, 3.0},        Ion{"LI", 3, 1.0},
    Ion{"LU", 71, 3.0},       Ion{"MN", 25, 2.0},        Ion{"MN3", 25, 3.0},
    Ion{"ND", 60, 3.0},       Ion{"NI", 28, 2.0},        Ion{"OS", 76, 3.0},
    Ion{"OS4", 76, 4.0},      Ion{"PB", 82, 2.0},        Ion{"PD", 46, 2.0},
    Ion{"PR", 59, 3.0},       Ion{"PT", 78, 2.0},        Ion{"PT4", 78, 4.0},
    Ion{"RB", 37, 1.0},       Ion{"RH", 45, 1.0},        Ion{"RH3", 45, 3.0},
    Ion{"RHF", 45, 2.0},      Ion{"RU", 44, 3.0},        Ion{"SB", 51, 3.0},
    Ion{"SM", 62, 3.0},       Ion{"SR", 38, 2.0},        Ion{"TB", 65, 3.0},
    Ion{"TH", 90, 4.0},       Ion{"TL", 81, 1.0},        Ion{"V", 23, 3.0},
    Ion{"W", 74, 6.0},        Ion{"Y1", 39, 2.0},        Ion{"YB", 70, 3.0},
    Ion{"YB2", 70, 2.0},      Ion{"YT3", 39, 3.0},       Ion{"ZCM", 96, 3.0},
    Ion{"ZR", 40, 4.0},       Ion{"ZTM", 89, 3.0},
};

auto context(const std::size_t molecule_index, const std::size_t atom_index,
             const std::string_view component, const std::string_view message) -> std::string {
    auto result = std::string{message} + " at record " + std::to_string(molecule_index);
    if (atom_index != std::string_view::npos) {
        result += ", atom " + std::to_string(atom_index);
    }
    if (!component.empty()) {
        result += " (component " + std::string{component} + ")";
    }
    return result;
}

} // namespace

auto fixed_ion_names(const bool common_only) -> std::vector<std::string> {
    std::vector<std::string> names;
    for (const auto& ion : ions) {
        if (!common_only || ion.common) {
            names.emplace_back(ion.component);
        }
    }
    if (!common_only) {
        std::ranges::sort(names);
    }
    return names;
}

auto resolve_fixed_ions(const std::span<const ImportedMoleculeRecord> records,
                        const std::span<const std::string> component_ids)
    -> calculation::FixedIons {
    calculation::FixedIons result;
    if (component_ids.empty()) {
        return result;
    }

    std::vector<const Ion*> selected;
    for (const auto& id : component_ids) {
        const auto ion = std::ranges::find(ions, id, &Ion::component);
        if (ion == ions.end()) {
            throw std::invalid_argument("unknown fixed-ion component ID: " + id);
        }
        if (std::ranges::find(selected, std::addressof(*ion)) == selected.end()) {
            selected.push_back(std::addressof(*ion));
        }
    }

    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const auto& record = records[molecule_index];
        if (record.molecule.atom_count() == 0) {
            continue;
        }
        if (!record.import_metadata || record.import_metadata->components.empty()) {
            throw std::invalid_argument(
                "fixed ions are selected by residue name, but " +
                core::detail::molecule_description(record.molecule, molecule_index) +
                " has no residues; read the structure from a PDB or mmCIF file");
        }
        if (record.import_metadata->atoms.size() != record.molecule.atom_count()) {
            throw std::invalid_argument(
                context(molecule_index, std::string_view::npos, {},
                        "fixed-ion source atom mapping does not match the molecule's atoms"));
        }
        const auto& atoms = record.import_metadata->atoms;
        const auto atom_count = record.molecule.atom_count();
        constexpr auto missing_component = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> component_for_atom(atom_count, missing_component);
        std::vector<bool> bonds(record.molecule.atom_count(), false);
        for (const auto& bond : record.molecule.bonds()) {
            bonds[bond.first_atom_index()] = true;
            bonds[bond.second_atom_index()] = true;
        }
        for (std::size_t component_index = 0;
             component_index < record.import_metadata->components.size(); ++component_index) {
            const auto& component = record.import_metadata->components[component_index];
            if (component.component_id.empty() || component.atom_indices.empty()) {
                throw std::invalid_argument(context(molecule_index, std::string_view::npos,
                                                    component.component_id,
                                                    "component metadata entry is empty"));
            }
            for (const auto atom_index : component.atom_indices) {
                if (atom_index >= atom_count ||
                    component_for_atom[atom_index] != missing_component) {
                    throw std::invalid_argument(
                        context(molecule_index, atom_index, component.component_id,
                                "component metadata has an invalid or repeated atom index"));
                }
                component_for_atom[atom_index] = component_index;
            }
        }
        if (std::ranges::find(component_for_atom, missing_component) != component_for_atom.end()) {
            throw std::invalid_argument(context(molecule_index, std::string_view::npos, {},
                                                "component metadata does not cover every atom"));
        }

        const auto first_source = result.sources.size();
        for (const auto& component : record.import_metadata->components) {
            const auto canonical_name = canonical_component_name(
                atoms[component.atom_indices.front()].structural_labels, component.component_id);
            const auto target = std::ranges::find(selected, canonical_name,
                                                  [](const Ion* ion) { return ion->component; });
            if (target == selected.end()) {
                continue;
            }
            if (component.atom_indices.size() != 1) {
                throw std::invalid_argument(
                    context(molecule_index, component.atom_indices.front(), component.component_id,
                            "selected fixed-ion component is not monatomic"));
            }
            const auto atom_index = component.atom_indices.front();
            const auto& atom = record.molecule.atom(atom_index);
            if (atom.atomic_number() != (*target)->atomic_number) {
                throw std::invalid_argument(
                    context(molecule_index, atom_index, component.component_id,
                            "selected fixed-ion component atom identity mismatch"));
            }
            if (bonds[atom_index]) {
                throw std::invalid_argument(
                    context(molecule_index, atom_index, component.component_id,
                            "selected fixed-ion component has an incident graph bond"));
            }
            result.sources.push_back({molecule_index, atom_index, (*target)->charge});
        }
        std::ranges::sort(result.sources.begin() + static_cast<std::ptrdiff_t>(first_source),
                          result.sources.end(), {}, &calculation::FixedAtomCharge::atom_index);
    }
    return result;
}

} // namespace chargefw::adapters::detail
