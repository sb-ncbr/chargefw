#include "support/test_calculation.h"

#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <snitch/snitch.hpp>
#include <string>
#include <utility>
#include <vector>

namespace adapters = chargefw::adapters;
namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace parameters = chargefw::parameters;
namespace sdf = chargefw::adapters::native::sdf_input;

using chargefw::test::calculate_application;

namespace {

[[nodiscard]] auto representative_fixture() -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" / "corpus" /
           "sdf" / "audit" / "representative.sdf";
}

[[nodiscard]] auto representative_eem_parameters() -> parameters::ParameterSet {
    auto entries = std::vector<parameters::AtomParameterEntry>{};
    for (const auto atomic_number : {1, 6, 7, 8, 15, 16, 17}) {
        entries.push_back(
            {.key = {.atomic_number = atomic_number,
                     .classification = parameters::AtomParameterClassificationKind::plain,
                     .type = "*"},
             .parameters = {{.name = "A", .value = static_cast<double>(atomic_number)},
                            {.name = "B", .value = 20.0}}});
    }

    return parameters::ParameterSet{parameters::ParameterSetMetadata{.id = "representative-eem",
                                                                     .method_id = "eem",
                                                                     .name = "Representative EEM"},
                                    parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
                                    parameters::AtomParameters{std::move(entries)}};
}

} // namespace

TEST_CASE("representative molecules preserve full execution under whole-radius reduction",
          "[calculation][reduced-execution]") {
    std::ifstream input{representative_fixture()};
    REQUIRE(input);
    auto reader = sdf::SdfReader{input, "representative.sdf"};
    const auto parameter_sets = std::vector{representative_eem_parameters()};
    auto record_count = std::size_t{0};

    while (const auto record = reader.next()) {
        CAPTURE(record->molecule.name());
        const auto molecules = core::MoleculeCollection{std::vector{record->molecule}};
        const auto full = calculate_application(
            {.molecules = molecules, .parameter_sets = parameter_sets, .method_id = "eem"});
        REQUIRE(full.calculated());
        REQUIRE(full.charges.has_value());
        const auto& full_charges = full.charges->assignment(0).charges;
        CHECK(std::abs(full_charges.total() - core::total_formal_charge(record->molecule)) <
              1.0e-10);

        for (const auto selection : {calculation::ExecutionSelectionKind::cutoff,
                                     calculation::ExecutionSelectionKind::cover}) {
            const auto reduced = calculate_application(
                {.molecules = molecules,
                 .parameter_sets = parameter_sets,
                 .method_id = "eem",
                 .execution_selection = calculation::ExecutionSelection{selection, 20.0}});
            REQUIRE(reduced.calculated());
            REQUIRE(reduced.charges.has_value());
            chargefw::test::assert_same_charges(reduced.charges->assignment(0).charges,
                                                full_charges, 1.0e-10);
        }

        ++record_count;
    }

    CHECK(record_count == 8);
}
