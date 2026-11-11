#include "calculation/calculation_request.h"
#include "methods/find_applicable_methods.h"
#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <array>
#include <cmath>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

using Complexity = methods::ComplexityTerm;

#ifndef CHARGEFW_TEST_PARAMETER_DIR
#error "CHARGEFW_TEST_PARAMETER_DIR must be defined"
#endif

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared_molecule{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{prepared_molecule, method_options,
                                          chargefw::core::total_formal_charge(molecule)};

    return method.calculate(input);
}

auto make_smpqeq_water_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-smpqeq", .method_id = "smpqeq", .name = "Test SMP/QEq water parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "first", .value = 1.0},
                                                    {.name = "second", .value = 10.0},
                                                    {.name = "third", .value = 1.0},
                                                    {.name = "fourth", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "first", .value = 2.0},
                                                    {.name = "second", .value = 10.0},
                                                    {.name = "third", .value = 1.0},
                                                    {.name = "fourth", .value = 0.1}}}}}};
}

auto make_ammonium() -> chargefw::core::Molecule {
    std::vector atoms{chargefw::core::Atom{7, 1, "N"}, chargefw::core::Atom{1, 0, "H1"},
                      chargefw::core::Atom{1, 0, "H2"}, chargefw::core::Atom{1, 0, "H3"},
                      chargefw::core::Atom{1, 0, "H4"}};
    std::vector bonds{chargefw::core::Bond{0, 1, chargefw::core::BondOrder::SINGLE},
                      chargefw::core::Bond{0, 2, chargefw::core::BondOrder::SINGLE},
                      chargefw::core::Bond{0, 3, chargefw::core::BondOrder::SINGLE},
                      chargefw::core::Bond{0, 4, chargefw::core::BondOrder::SINGLE}};
    std::vector positions{chargefw::core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                          chargefw::core::Position{.x = 1.0, .y = 1.0, .z = 1.0},
                          chargefw::core::Position{.x = -1.0, .y = -1.0, .z = 1.0},
                          chargefw::core::Position{.x = -1.0, .y = 1.0, .z = -1.0},
                          chargefw::core::Position{.x = 1.0, .y = -1.0, .z = -1.0}};
    std::vector conformers{chargefw::core::Conformer{std::move(positions), "model-1"}};

    return chargefw::core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                                    "ammonium"};
}

auto make_formally_charged_neutral_pair() -> chargefw::core::Molecule {
    return chargefw::core::Molecule{
        std::vector{chargefw::core::Atom{7, 1, "N"}, chargefw::core::Atom{8, -1, "O"}},
        std::vector{chargefw::core::Bond{0, 1, chargefw::core::BondOrder::SINGLE}},
        {},
        "formal-charge-pair"};
}

auto make_component_test_molecule(std::vector<chargefw::core::Atom> atoms,
                                  std::vector<chargefw::core::Bond> bonds, std::string name)
    -> chargefw::core::Molecule {
    std::vector<chargefw::core::Position> positions;
    positions.reserve(atoms.size());
    for (std::size_t atom_index = 0; atom_index < atoms.size(); ++atom_index) {
        positions.push_back(
            chargefw::core::Position{.x = static_cast<double>(atom_index), .y = 0.0, .z = 0.0});
    }

    std::vector conformers{chargefw::core::Conformer{std::move(positions), std::string{"model-1"}}};
    return chargefw::core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                                    std::move(name)};
}

constexpr std::array component_neutral_method_ids{
    std::string_view{"charge2"}, std::string_view{"delre"}, std::string_view{"denr"},
    std::string_view{"gdac"},    std::string_view{"kcm"},   std::string_view{"mgc"},
    std::string_view{"mpeoe"},   std::string_view{"peoe"},  std::string_view{"sqe"}};

// Checks prerequisites of each method on the molecule. A rejection must be a single
// unsupported-molecule issue whose message contains the given context.
auto check_neutrality_prerequisites(const chargefw::core::Molecule& molecule,
                                    const std::span<const std::string_view> method_ids,
                                    const std::optional<std::string_view> rejection_context)
    -> void {
    const features::PreparedMolecule prepared{molecule};
    for (const auto method_id : method_ids) {
        CAPTURE(method_id);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);

        const auto options = methods::make_default_options(method->option_schema());
        const auto result = method->check_method_prerequisites(
            {.prepared_molecule = prepared, .method_options = options});
        CHECK(static_cast<bool>(result) == !rejection_context.has_value());

        if (rejection_context.has_value()) {
            REQUIRE(result.issues.size() == 1);
            CHECK(result.issues[0].kind == methods::PrerequisiteIssueKind::unsupported_molecule);
            CHECK(result.issues[0].message.contains(*rejection_context));
        }
    }
}

} // namespace

TEST_CASE("method complexity terms have Big-O notation", "[methods]") {
    CHECK(methods::complexity_notation(Complexity::constant) == "O(1)");
    CHECK(methods::complexity_notation(Complexity::atoms_cubed) == "O(n^3)");
    CHECK(methods::complexity_notation(Complexity::atoms_squared_plus_bonds_squared) ==
          "O(n^2 + m^2)");
    CHECK(methods::complexity_notation(Complexity::atoms_plus_bonds_cubed) == "O((n + m)^3)");
}

TEST_CASE("built-in resource estimates grow with their declared complexity",
          "[methods][builtin-methods]") {
    const auto degree = [](const Complexity term) -> double {
        switch (term) {
        case Complexity::constant:
            return 0.0;
        case Complexity::atoms:
        case Complexity::bonds:
        case Complexity::atoms_plus_bonds:
            return 1.0;
        case Complexity::atoms_squared:
        case Complexity::bonds_squared:
        case Complexity::atoms_plus_bonds_squared:
        case Complexity::atoms_squared_plus_bonds_squared:
            return 2.0;
        case Complexity::atoms_cubed:
        case Complexity::bonds_cubed:
        case Complexity::atoms_plus_bonds_cubed:
            return 3.0;
        }
        return -1.0;
    };
    const auto growth = [](const double small, const double large) -> double {
        return small == 0.0 && large == 0.0 ? 0.0 : std::log2(large / small);
    };
    // Joint scaling keeps the mean degree fixed and makes lower-order terms negligible.
    constexpr auto size = 1.0e5;

    for (const auto& method : methods::method_registry().methods()) {
        CAPTURE(method->id());
        const auto options = methods::make_default_options(method->option_schema());
        const auto small = method->estimate({.atoms = size, .bonds = size}, options);
        const auto large = method->estimate({.atoms = 2.0 * size, .bonds = 2.0 * size}, options);
        const auto resources = method->requirements().resources;
        CHECK(std::abs(growth(small.operations, large.operations) - degree(resources.time)) < 0.05);
        CHECK(std::abs(growth(small.memory_bytes, large.memory_bytes) - degree(resources.memory)) <
              0.05);
    }
}

TEST_CASE("every built-in completes its declared full workflow", "[methods][builtin-methods]") {
    const auto molecule = chargefw::test::make_two_conformer_water();
    const auto collection = chargefw::core::MoleculeCollection{std::vector{molecule}, "water"};
    const auto prepared = features::PreparedMoleculeCollection{collection};
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    const std::array smpqeq_parameter_sets{make_smpqeq_water_parameters()};

    REQUIRE_FALSE(methods::method_registry().methods().empty());
    for (const auto& method : methods::method_registry().methods()) {
        CAPTURE(method->id());

        const std::array<const methods::Method*, 1> candidates{method.get()};
        const auto candidate_parameter_sets =
            method->id() == "smpqeq"
                ? std::span<const parameters::ParameterSet>{smpqeq_parameter_sets}
                : std::span<const parameters::ParameterSet>{parameter_sets};
        const auto applicability =
            methods::find_applicable_methods({.molecules = prepared,
                                              .methods = candidates,
                                              .parameter_sets = candidate_parameter_sets});
        REQUIRE_FALSE(applicability.applicable.empty());
        const auto& selected = applicability.applicable.front();

        const auto result =
            chargefw::calculation::calculate({.molecules = prepared, .selected = selected});
        CHECK(result.charges.method_id() == method->id());
        CHECK(result.charges.parameter_set_id().has_value() == method->requires_parameters());
        CHECK(result.charges.size() ==
              (method->requirements().coordinates ? molecule.conformer_count() : 1));
        if (method->requirements().coordinates) {
            REQUIRE(result.charges.size() == 2);
            CHECK(std::abs(result.charges.assignment(0).charges[0] -
                           result.charges.assignment(1).charges[0]) > 1.0e-4);
        }

        for (const auto& assignment : result.charges.assignments()) {
            CHECK(assignment.target.molecule_index == 0);
            CHECK(assignment.charges.size() == molecule.atom_count());
            for (const auto charge : assignment.charges.values()) {
                CHECK(std::isfinite(charge));
            }
            CHECK(std::abs(assignment.charges.total()) < 1.0e-8);
        }

        if (!method->requirements().coordinates) {
            const auto graph_water = chargefw::test::make_water_graph();
            const auto graph_collection =
                chargefw::core::MoleculeCollection{std::vector{graph_water}, "water"};
            const auto graph_prepared = features::PreparedMoleculeCollection{graph_collection};
            const auto graph_applicability =
                methods::find_applicable_methods({.molecules = graph_prepared,
                                                  .methods = candidates,
                                                  .parameter_sets = candidate_parameter_sets});
            REQUIRE_FALSE(graph_applicability.applicable.empty());
            const auto& graph_selected = graph_applicability.applicable.front();

            const auto graph_result = chargefw::calculation::calculate(
                {.molecules = graph_prepared, .selected = graph_selected});
            REQUIRE(graph_result.charges.size() == 1);
            chargefw::test::assert_same_charges(graph_result.charges.assignment(0).charges,
                                                result.charges.assignment(0).charges, 1.0e-10);
        }

        if (method->id() != "formal" && method->id() != "dummy") {
            chargefw::test::assert_neutral_water_charges(result.charges.assignment(0).charges,
                                                         1.0e-10);
        }

        const std::vector<parameters::ParameterSet> selected_parameter_sets =
            selected.parameter_set == nullptr
                ? std::vector<parameters::ParameterSet>{}
                : std::vector<parameters::ParameterSet>{*selected.parameter_set};
        chargefw::test::assert_water_charges_labeling_invariant(
            method->id(), result.charges.assignment(0).charges, selected_parameter_sets);
    }
}

TEST_CASE("neutral-only methods reject ammonium with bundled parameters",
          "[methods][builtin-methods]") {
    const auto ammonium = make_ammonium();
    const auto collection = chargefw::core::MoleculeCollection{std::vector{ammonium}, "ammonium"};
    const auto prepared = features::PreparedMoleculeCollection{collection};
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    constexpr std::array method_ids{std::string_view{"peoe"}, std::string_view{"mpeoe"},
                                    std::string_view{"delre"}, std::string_view{"sqe"}};

    for (const auto method_id : method_ids) {
        CAPTURE(method_id);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);

        const std::array candidates{method};
        const auto applicability = methods::find_applicable_methods(
            {.molecules = prepared, .methods = candidates, .parameter_sets = parameter_sets});

        CHECK(applicability.applicable.empty());
        REQUIRE(applicability.rejected.size() == 1);
        CHECK_FALSE(applicability.rejected[0].parameter_set_index.has_value());
        REQUIRE(applicability.rejected[0].issues.size() == 1);
        const auto& issue = applicability.rejected[0].issues[0];
        CHECK(issue.kind == methods::PrerequisiteIssueKind::unsupported_molecule);
        CHECK(issue.message.contains("zero total formal charge"));
    }
}

TEST_CASE("zero-initialized methods reject a cancelling disconnected ion pair",
          "[methods][builtin-methods]") {
    const auto molecule = make_component_test_molecule(
        {chargefw::core::Atom{7, 1, "N"}, chargefw::core::Atom{1, 0, "H1"},
         chargefw::core::Atom{1, 0, "H2"}, chargefw::core::Atom{1, 0, "H3"},
         chargefw::core::Atom{1, 0, "H4"}, chargefw::core::Atom{8, -1, "O"},
         chargefw::core::Atom{1, 0, "H5"}},
        {chargefw::core::Bond{0, 1, chargefw::core::BondOrder::SINGLE},
         chargefw::core::Bond{0, 2, chargefw::core::BondOrder::SINGLE},
         chargefw::core::Bond{0, 3, chargefw::core::BondOrder::SINGLE},
         chargefw::core::Bond{0, 4, chargefw::core::BondOrder::SINGLE},
         chargefw::core::Bond{5, 6, chargefw::core::BondOrder::SINGLE}},
        "ammonium-hydroxide-pair");

    CHECK(chargefw::core::total_formal_charge(molecule) == 0.0);
    check_neutrality_prerequisites(molecule, component_neutral_method_ids,
                                   "every bonded part has zero total formal charge");
}

TEST_CASE("zero-initialized methods accept neutral disconnected components",
          "[methods][builtin-methods]") {
    const auto molecule = make_component_test_molecule(
        {chargefw::core::Atom{8, 0, "O1"}, chargefw::core::Atom{1, 0, "H1"},
         chargefw::core::Atom{8, 0, "O2"}, chargefw::core::Atom{1, 0, "H2"}},
        {chargefw::core::Bond{0, 1, chargefw::core::BondOrder::SINGLE},
         chargefw::core::Bond{2, 3, chargefw::core::BondOrder::SINGLE}},
        "neutral-components");

    check_neutrality_prerequisites(molecule, component_neutral_method_ids, std::nullopt);
}

TEST_CASE("zero-initialized methods accept a connected net-neutral zwitterion",
          "[methods][builtin-methods]") {
    const auto molecule = make_component_test_molecule(
        {chargefw::core::Atom{7, 1, "N"}, chargefw::core::Atom{6, 0, "C"},
         chargefw::core::Atom{8, -1, "O"}},
        {chargefw::core::Bond{0, 1, chargefw::core::BondOrder::SINGLE},
         chargefw::core::Bond{1, 2, chargefw::core::BondOrder::SINGLE}},
        "zwitterion");

    check_neutrality_prerequisites(molecule, component_neutral_method_ids, std::nullopt);
}

TEST_CASE("neutral-only methods reject a cation", "[methods][builtin-methods]") {
    const auto cation =
        make_component_test_molecule({chargefw::core::Atom{1, 1, "H"}}, {}, "cation");

    check_neutrality_prerequisites(cation, component_neutral_method_ids,
                                   "every bonded part has zero total formal charge");
    check_neutrality_prerequisites(cation, std::array{std::string_view{"veem"}},
                                   "neutral molecules");
}

TEST_CASE("built-in iteration options reach the charge calculation", "[methods][builtin-methods]") {
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    auto denr_options = methods::MethodOptions{};
    denr_options.set("step", 0.1);
    denr_options.set("iters", 25);
    auto charge2_options = methods::MethodOptions{};
    charge2_options.set("iters", 1);
    // Charge2 polarizability feedback needs a chain longer than water to respond to iterations.
    const auto hco_chain =
        chargefw::core::Molecule{{chargefw::core::Atom{1, 0, "H"}, chargefw::core::Atom{6, 0, "C"},
                                  chargefw::core::Atom{8, 0, "O"}},
                                 {chargefw::core::Bond{0, 1}, chargefw::core::Bond{1, 2}},
                                 {},
                                 "HCO chain"};
    struct OptionCase {
        std::string_view method_id;
        chargefw::core::Molecule molecule;
        methods::MethodOptions options;
    };
    const auto cases =
        std::array{OptionCase{"denr", chargefw::test::make_water_graph(), denr_options},
                   OptionCase{"charge2", hco_chain, charge2_options}};

    for (const auto& [method_id, molecule, options] : cases) {
        CAPTURE(method_id);
        const auto collection = chargefw::core::MoleculeCollection{std::vector{molecule}};
        const auto prepared = features::PreparedMoleculeCollection{collection};
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        const std::array candidates{method};
        const auto calculate_with = [&](const methods::MethodOptions* configured) {
            auto method_options = std::unordered_map<std::string, methods::MethodOptions>{};
            if (configured != nullptr) {
                method_options.emplace(std::string{method_id}, *configured);
            }
            const auto applicability =
                methods::find_applicable_methods({.molecules = prepared,
                                                  .methods = candidates,
                                                  .parameter_sets = parameter_sets,
                                                  .method_options = std::move(method_options)});
            REQUIRE_FALSE(applicability.applicable.empty());
            return chargefw::calculation::calculate(
                       {.molecules = prepared, .selected = applicability.applicable.front()})
                .charges.assignment(0)
                .charges;
        };

        const auto default_charges = calculate_with(nullptr);
        const auto configured_charges = calculate_with(&options);
        CHECK(std::abs(default_charges[0] - configured_charges[0]) > 1.0e-6);
        CHECK(std::abs(configured_charges.total()) < 1.0e-10);
    }
}

TEST_CASE("PEOE methods initialize from formal charges when requested",
          "[methods][builtin-methods]") {
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    constexpr std::array method_ids{std::string_view{"peoe"}, std::string_view{"mpeoe"}};

    for (const auto method_id : method_ids) {
        CAPTURE(method_id);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        const std::array candidates{method};
        auto formal_options = methods::MethodOptions{};
        formal_options.set("initial_charges", std::string{"formal"});

        const auto ammonium = make_ammonium();
        const auto ammonium_collection =
            chargefw::core::MoleculeCollection{std::vector{ammonium}, "ammonium"};
        const auto prepared_ammonium = features::PreparedMoleculeCollection{ammonium_collection};
        const auto ammonium_applicability = methods::find_applicable_methods(
            {.molecules = prepared_ammonium,
             .methods = candidates,
             .parameter_sets = parameter_sets,
             .method_options = {{std::string{method_id}, formal_options}}});

        CHECK(ammonium_applicability.rejected.empty());
        REQUIRE(ammonium_applicability.applicable.size() == 1);
        const auto ammonium_result = chargefw::calculation::calculate(
            {.molecules = prepared_ammonium,
             .selected = ammonium_applicability.applicable.front()});
        REQUIRE(ammonium_result.charges.size() == 1);
        CHECK(std::abs(ammonium_result.charges.assignment(0).charges.total() - 1.0) < 1.0e-10);

        const auto neutral_pair = make_formally_charged_neutral_pair();
        const auto pair_collection =
            chargefw::core::MoleculeCollection{std::vector{neutral_pair}, "formal-charge-pair"};
        const auto prepared_pair = features::PreparedMoleculeCollection{pair_collection};
        const auto zero_applicability = methods::find_applicable_methods(
            {.molecules = prepared_pair, .methods = candidates, .parameter_sets = parameter_sets});
        const auto formal_applicability = methods::find_applicable_methods(
            {.molecules = prepared_pair,
             .methods = candidates,
             .parameter_sets = parameter_sets,
             .method_options = {{std::string{method_id}, formal_options}}});

        REQUIRE(zero_applicability.applicable.size() == 1);
        REQUIRE(formal_applicability.applicable.size() == 1);
        const auto zero_result = chargefw::calculation::calculate(
            {.molecules = prepared_pair, .selected = zero_applicability.applicable.front()});
        const auto formal_result = chargefw::calculation::calculate(
            {.molecules = prepared_pair, .selected = formal_applicability.applicable.front()});
        const auto& zero_charges = zero_result.charges.assignment(0).charges;
        const auto& formal_charges = formal_result.charges.assignment(0).charges;

        CHECK(std::abs(zero_charges.total()) < 1.0e-10);
        CHECK(std::abs(formal_charges.total()) < 1.0e-10);
        CHECK(std::abs(zero_charges[0] - formal_charges[0]) > 1.0e-6);
    }
}

TEST_CASE("dummy assigns exact zero charges", "[methods][builtin-methods]") {
    const auto* dummy = methods::method_registry().find("dummy");
    REQUIRE(dummy != nullptr);

    const auto water = chargefw::test::make_water_graph();
    const auto dummy_charges = calculate(*dummy, water);
    CHECK(dummy_charges.size() == water.atom_count());

    for (const auto charge : dummy_charges.values()) {
        CHECK(charge == 0.0);
    }
}

TEST_CASE("formal copies atomic formal charges", "[methods][builtin-methods]") {
    const auto* formal = methods::method_registry().find("formal");
    REQUIRE(formal != nullptr);

    const auto charged_pair = chargefw::test::make_formally_charged_pair();
    const auto formal_charges = calculate(*formal, charged_pair);
    REQUIRE(formal_charges.size() == charged_pair.atom_count());
    CHECK(formal_charges[0] == 1.0);
    CHECK(formal_charges[1] == -1.0);
    CHECK(formal_charges.total() == 0.0);
}
