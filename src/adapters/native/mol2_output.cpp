#include <chargefw/adapters/native/mol2_output.h>

#include <chargefw/core/periodic_table.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <memory>
#include <ostream>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace chargefw::adapters::native::mol2_output {
namespace {

[[nodiscard]] auto safe_line(std::string value) -> std::string {
    std::ranges::replace_if(
        value, [](const char character) { return character == '\r' || character == '\n'; }, '_');
    return value;
}

[[nodiscard]] auto record_name(const ImportedMoleculeRecord& record,
                               const std::size_t molecule_index, const std::size_t conformer_index)
    -> std::string {
    auto result = record.identity.record_id.empty() ? std::string{record.molecule.name()}
                                                    : record.identity.record_id.display_string();
    if (result.empty()) {
        result = "molecule_" + std::to_string(molecule_index + 1);
    }
    result = safe_line(std::move(result));
    if (record.molecule.conformer_count() > 1) {
        result += "_conformer_" + std::to_string(conformer_index + 1);
    }
    return result;
}

[[nodiscard]] auto atom_name(const core::Atom& atom, const std::size_t atom_index) -> std::string {
    const auto name = atom.name();
    if (!name.empty() && std::ranges::none_of(name, [](const unsigned char character) {
            return std::isspace(character) != 0;
        })) {
        return std::string{name};
    }
    return std::string{core::element_symbol(atom.atomic_number())} + std::to_string(atom_index + 1);
}

[[nodiscard]] auto assignment_for(const charges::ChargeSet& charge_set,
                                  const std::size_t molecule_index,
                                  const std::size_t conformer_index)
    -> const charges::ChargeAssignment& {
    const auto found = std::ranges::find_if(
        charge_set.assignments(), [molecule_index, conformer_index](const auto& assignment) {
            return assignment.target.molecule_index == molecule_index &&
                   (!assignment.target.conformer_index.has_value() ||
                    *assignment.target.conformer_index == conformer_index);
        });
    if (found == charge_set.assignments().end()) {
        throw std::invalid_argument{"MOL2 output has no charge assignment for molecule " +
                                    std::to_string(molecule_index + 1) + ", conformer " +
                                    std::to_string(conformer_index + 1)};
    }
    return *found;
}

auto validate_result(const ChargeCalculationResult& result) -> void {
    if (!result.execution().calculated() || !result.execution().charges.has_value()) {
        throw std::invalid_argument{"MOL2 output requires a successful calculation"};
    }
    if (result.inputs().empty()) {
        throw std::invalid_argument{"MOL2 output requires at least one molecule record"};
    }

    const auto& charge_set = *result.execution().charges;
    for (std::size_t molecule_index = 0; molecule_index < result.inputs().size();
         ++molecule_index) {
        const auto& molecule = result.inputs()[molecule_index].molecule;
        if (molecule.atom_count() == 0) {
            throw std::invalid_argument{"MOL2 output requires at least one atom per molecule"};
        }
        if (molecule.conformer_count() == 0) {
            throw std::invalid_argument{"MOL2 output requires coordinates for molecule " +
                                        std::to_string(molecule_index + 1)};
        }
        for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
             ++conformer_index) {
            const auto& conformer = molecule.conformer(conformer_index);
            const auto& assignment = assignment_for(charge_set, molecule_index, conformer_index);
            if (assignment.charges.size() != molecule.atom_count()) {
                throw std::invalid_argument{
                    "MOL2 charge assignment size does not match molecule atom count"};
            }
            for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
                const auto& position = conformer[atom_index];
                if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                    !std::isfinite(position.z)) {
                    throw std::invalid_argument{"MOL2 coordinates must be finite"};
                }
                if (!std::isfinite(assignment.charges[atom_index])) {
                    throw std::invalid_argument{"MOL2 charges must be finite"};
                }
            }
        }
    }
}

auto write_record(std::ostream& output, const ImportedMoleculeRecord& record,
                  const std::size_t molecule_index, const std::size_t conformer_index,
                  const charges::ChargeAssignment& assignment, const charges::ChargeSet& charge_set,
                  const std::string_view generator_version) -> void {
    const auto& molecule = record.molecule;
    const auto& conformer = molecule.conformer(conformer_index);
    const auto name = record_name(record, molecule_index, conformer_index);
    std::println(output, "@<TRIPOS>MOLECULE");
    std::println(output, "{}", name);
    std::println(output, "{} {} 1 0 0", molecule.atom_count(), molecule.bond_count());
    std::println(output, "SMALL");
    std::println(output, "USER_CHARGES");
    std::println(output, "****");
    std::print(output, "Generated by ChargeFW");
    if (!generator_version.empty()) {
        std::print(output, " {}", safe_line(std::string{generator_version}));
    }
    std::print(output, "; method={}", safe_line(std::string{charge_set.method_id()}));
    if (const auto parameter_set_id = charge_set.parameter_set_id(); parameter_set_id.has_value()) {
        std::print(output, "; parameter_set={}", safe_line(std::string{*parameter_set_id}));
    }
    std::println(output);
    std::println(output, "@<TRIPOS>ATOM");
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);
        const auto& position = conformer[atom_index];
        std::println(output, "{} {} {} {} {} {} 1 UNL {:.4f}", atom_index + 1,
                     atom_name(atom, atom_index), position.x, position.y, position.z,
                     core::element_symbol(atom.atomic_number()), assignment.charges[atom_index]);
    }

    std::println(output, "@<TRIPOS>BOND");
    for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        std::println(output, "{} {} {} {}", bond_index + 1, bond.first_atom_index() + 1,
                     bond.second_atom_index() + 1, core::to_string(bond.order()));
    }
    std::println(output, "@<TRIPOS>SUBSTRUCTURE");
    std::println(output, "1 UNL 1");
}

} // namespace

Mol2Writer::Mol2Writer(std::ostream& output) : output_{std::addressof(output)} {}

auto Mol2Writer::write(const ChargeCalculationResult& result,
                       const std::string_view generator_version) const -> void {
    validate_result(result);
    const auto& charge_set = *result.execution().charges;
    bool first = true;
    for (std::size_t molecule_index = 0; molecule_index < result.inputs().size();
         ++molecule_index) {
        const auto& record = result.inputs()[molecule_index];
        for (std::size_t conformer_index = 0; conformer_index < record.molecule.conformer_count();
             ++conformer_index) {
            if (!first) {
                std::println(*output_);
            }
            first = false;
            write_record(*output_, record, molecule_index, conformer_index,
                         assignment_for(charge_set, molecule_index, conformer_index), charge_set,
                         generator_version);
        }
    }
    if (!*output_) {
        throw std::runtime_error{"failed to write MOL2 output"};
    }
}

} // namespace chargefw::adapters::native::mol2_output
