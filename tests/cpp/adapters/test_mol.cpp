#include <chargefw/adapters/native/mol2_input.h>
#include <chargefw/adapters/native/mol_input.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/core/bond.h>
#include <exception>
#include <filesystem>
#include <fstream>
#include <istream>
#include <locale>
#include <optional>
#include <snitch/snitch.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace adapters = chargefw::adapters;
namespace mol = chargefw::adapters::native::mol_input;
namespace mol2 = chargefw::adapters::native::mol2_input;
namespace sdf = chargefw::adapters::native::sdf_input;

static_assert(!std::is_copy_constructible_v<mol::MolReader> &&
              !std::is_copy_assignable_v<mol::MolReader>);
static_assert(!std::is_copy_constructible_v<sdf::SdfReader> &&
              !std::is_copy_assignable_v<sdf::SdfReader>);
static_assert(!std::is_copy_constructible_v<mol2::Mol2Reader> &&
              !std::is_copy_assignable_v<mol2::Mol2Reader>);
static_assert(std::is_move_constructible_v<mol::MolReader> &&
              std::is_move_constructible_v<sdf::SdfReader> &&
              std::is_move_constructible_v<mol2::Mol2Reader>);

namespace {

[[nodiscard]] auto fixture(const std::string_view name) -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" / name;
}

[[nodiscard]] auto v3000_with_coordinates(const std::string_view coordinates) -> std::string {
    return std::string{"coordinates\nchargefw\n\n"
                       "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                       "M  V30 BEGIN CTAB\n"
                       "M  V30 COUNTS 1 0 0 0 0\n"
                       "M  V30 BEGIN ATOM\n"
                       "M  V30 1 C "} +
           std::string{coordinates} +
           " 0\n"
           "M  V30 END ATOM\n"
           "M  V30 BEGIN BOND\n"
           "M  V30 END BOND\n"
           "M  V30 END CTAB\n"
           "M  END\n";
}

[[nodiscard]] auto with_crlf(const std::string_view source) -> std::string {
    auto result = std::string{};
    result.reserve(source.size());
    for (const auto character : source) {
        if (character == '\n') {
            result.push_back('\r');
        }
        result.push_back(character);
    }
    return result;
}

class DecimalCommaNumpunct final : public std::numpunct<char> {
  protected:
    [[nodiscard]] auto do_decimal_point() const -> char override {
        return ',';
    }
};

class ScopedGlobalLocale final {
  public:
    explicit ScopedGlobalLocale(const std::locale& replacement) {
        std::locale::global(replacement);
    }

    ScopedGlobalLocale(const ScopedGlobalLocale&) = delete;
    auto operator=(const ScopedGlobalLocale&) -> ScopedGlobalLocale& = delete;
    ScopedGlobalLocale(ScopedGlobalLocale&&) = delete;
    auto operator=(ScopedGlobalLocale&&) -> ScopedGlobalLocale& = delete;

    ~ScopedGlobalLocale() {
        std::locale::global(previous_);
    }

  private:
    std::locale previous_{};
};

} // namespace

TEST_CASE("native numeric parsing is complete and locale independent", "[adapters][native]") {
    {
        std::istringstream input{v3000_with_coordinates("+1.25 -2.5e1 3.125E-2")};
        const auto record = mol::parse_mol(input, {});
        const auto& position = record.molecule.conformer(0)[0];
        CHECK(position.x == 1.25);
        CHECK(position.y == -25.0);
        CHECK(position.z == 0.03125);
    }

    for (const auto* const coordinates : {"1.25suffix 0 0", "++1.25 0 0"}) {
        std::istringstream input{v3000_with_coordinates(coordinates)};

        CHECK_THROWS_AS(mol::parse_mol(input, {}), std::runtime_error);
    }

    {
        const auto comma_locale = std::locale{std::locale::classic(), new DecimalCommaNumpunct{}};
        const auto locale_guard = ScopedGlobalLocale{comma_locale};
        std::istringstream input{v3000_with_coordinates("1.25 2.5e1 -3.125E-2")};
        const auto record = mol::parse_mol(input, {});
        const auto& position = record.molecule.conformer(0)[0];
        CHECK(position.x == 1.25);
        CHECK(position.y == 25.0);
        CHECK(position.z == -0.03125);
    }
}

TEST_CASE("native molecular mappings retain original atom tokens", "[adapters][mol]") {
    auto mol_record = [] {
        std::istringstream input{"mapping\nchargefw\n\n"
                                 "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                                 "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 0 0 0 0\nM  V30 BEGIN ATOM\n"
                                 "M  V30 001 C 0 0 0 0\nM  V30 10 O 1 0 0 0\nM  V30 END ATOM\n"
                                 "M  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"};
        auto reader = mol::MolReader{input, "mapping.mol"};
        return *reader.next();
    }();

    REQUIRE(mol_record.import_metadata.has_value());
    const auto& mol_mapping = *mol_record.import_metadata;
    CHECK(mol_mapping.format == adapters::MolecularSourceFormat::mol);
    CHECK(mol_mapping.atoms[0].id == "001");
    CHECK(mol_mapping.atoms[1].id == "10");
    REQUIRE(mol_mapping.conformers.size() == 1);
    CHECK(mol_mapping.conformers[0].sites == mol_mapping.atoms);
    CHECK(mol_mapping.source_connectivity == adapters::SourceConnectivity::explicitly_empty);

    std::istringstream mol2_input{
        "@<TRIPOS>MOLECULE\nmapping\n2 0 0 0 0\nSMALL\nNO_CHARGES\n"
        "@<TRIPOS>ATOM\n001 C1 0 0 0 C.3\n10 O1 1 0 0 O.2\n@<TRIPOS>BOND\n"};
    auto mol2_reader = mol2::Mol2Reader{mol2_input, "mapping.mol2"};
    const auto mol2_record = mol2_reader.next();
    REQUIRE(mol2_record.has_value());
    REQUIRE(mol2_record->import_metadata.has_value());
    CHECK(mol2_record->import_metadata->format == adapters::MolecularSourceFormat::mol2);
    CHECK(mol2_record->import_metadata->atoms[0].id == "001");
    CHECK(mol2_record->import_metadata->atoms[1].id == "10");
}

TEST_CASE("native MOL, SDF, and MOL2 input accepts CRLF", "[adapters][native]") {
    constexpr auto v2000 = "minimal\nchargefw\n\n"
                           "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                           "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                           "M  END\n";

    {
        auto input = std::istringstream{with_crlf(v2000)};
        const auto record = mol::parse_mol(input, {});
        CHECK(record.identity.record_id == "minimal");
        CHECK(record.molecule.atom_count() == 1);
    }

    {
        auto input = std::istringstream{with_crlf(v3000_with_coordinates("0 0 0"))};
        const auto record = mol::parse_mol(input, {});
        CHECK(record.identity.record_id == "coordinates");
        CHECK(record.molecule.atom_count() == 1);
    }

    {
        auto input =
            std::istringstream{with_crlf(std::string{v2000} + "$$$$\n" + v2000 + "$$$$\n")};
        auto reader = sdf::SdfReader{input, "records.sdf"};
        REQUIRE(reader.next().has_value());
        REQUIRE(reader.next().has_value());
        CHECK_FALSE(reader.next().has_value());
    }

    {
        constexpr auto untitled =
            "\nchargefw\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n";
        auto input = std::istringstream{
            with_crlf(std::string{v2000} + "$$$$\n" + untitled + "$$$$\n\n  \n")};
        input.exceptions(std::ios::failbit | std::ios::badbit);
        auto reader = sdf::SdfReader{input, "trailing_blank_lines.sdf"};
        REQUIRE(reader.next().has_value());
        const auto second = reader.next();
        REQUIRE(second.has_value());
        CHECK(second->molecule.atom_count() == 1);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        constexpr auto source = "@<TRIPOS>MOLECULE\nminimal\n2 1 0 0 0\nSMALL\nNO_CHARGES\n\n"
                                "@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 H1 1 0 0 H\n"
                                "@<TRIPOS>BOND\n1 1 2 1\n";
        auto lf_input = std::istringstream{source};
        auto crlf_input = std::istringstream{with_crlf(source)};
        auto lf_reader = mol2::Mol2Reader{lf_input, "lf.mol2"};
        auto crlf_reader = mol2::Mol2Reader{crlf_input, "crlf.mol2"};
        const auto lf_record = lf_reader.next();
        const auto crlf_record = crlf_reader.next();
        REQUIRE(lf_record.has_value());
        REQUIRE(crlf_record.has_value());
        CHECK(crlf_record->identity.record_id == lf_record->identity.record_id);
        CHECK(crlf_record->molecule.atom_count() == lf_record->molecule.atom_count());
        CHECK(crlf_record->molecule.bond_count() == lf_record->molecule.bond_count());
        CHECK(crlf_record->molecule.atom(0).name() == lf_record->molecule.atom(0).name());
        CHECK_FALSE(crlf_reader.next().has_value());
    }
}

TEST_CASE("native MOL input preserves V2000 and V3000 formal charges and bonds",
          "[adapters][native]") {
    {
        std::ifstream input{fixture("synthetic/mol/v2000/charged_atoms.mol")};
        auto reader = mol::MolReader{input, "charged_v2000.mol"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& record = *result;
        CHECK(record.molecule.atom_count() == 2);
        CHECK(record.molecule.bond_count() == 1);
        CHECK(record.molecule.conformer_count() == 1);
        CHECK(record.molecule.atom(0).formal_charge() == 1);
        CHECK(record.molecule.atom(1).formal_charge() == -1);
        // Molfiles have no atom names; the element symbol is not one.
        CHECK(record.molecule.atom(0).name().empty());
        CHECK(record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(record.identity.record_id == "charged-v2000");
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/mol/v3000/charged_atoms.mol")};
        const auto result = mol::parse_mol(
            input, {.source = "charged_v3000.mol", .record_index = 0, .record_id = {}});
        CHECK(result.molecule.atom_count() == 2);
        CHECK(result.molecule.atom(0).formal_charge() == 1);
        CHECK(result.molecule.atom(1).formal_charge() == -1);
        CHECK(result.molecule.atom(0).name().empty());
        CHECK(result.molecule.bond(0).first_atom_index() == 0);
        CHECK(result.molecule.bond(0).second_atom_index() == 1);
        CHECK(result.diagnostics.empty());
    }
}

TEST_CASE("native MOL and SDF input report ignored properties and bond orders",
          "[adapters][native]") {
    {
        std::ifstream input{fixture("synthetic/mol/v2000/ignored_properties.mol")};
        const auto result = mol::parse_mol(input, {});
        CHECK(result.molecule.atom_count() == 2);
        CHECK(result.molecule.bond_count() == 1);
        CHECK(result.molecule.atom(0).formal_charge() == 1);
        CHECK(result.molecule.atom(1).formal_charge() == -1);
        REQUIRE(result.diagnostics.size() == 2);
        CHECK(result.diagnostics[0].code == "v2000_property_ignored");
        CHECK(result.diagnostics[0].line == 8);
        CHECK(result.diagnostics[1].code == "v2000_property_ignored");
        CHECK(result.diagnostics[1].line == 9);
    }

    {
        std::ifstream input{fixture("synthetic/mol/v3000/bond_attributes_and_ignored_orders.mol")};
        const auto result = mol::parse_mol(input, {});
        CHECK(result.molecule.atom_count() == 4);
        REQUIRE(result.molecule.bond_count() == 2);
        CHECK(result.molecule.bond(0).first_atom_index() == 0);
        CHECK(result.molecule.bond(0).second_atom_index() == 1);
        CHECK(result.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(result.molecule.bond(1).first_atom_index() == 1);
        CHECK(result.molecule.bond(1).second_atom_index() == 3);
        CHECK(result.molecule.bond(1).order() == chargefw::core::BondOrder::DOUBLE);
        REQUIRE(result.diagnostics.size() == 2);
        CHECK(result.diagnostics[0].code == "v3000_bond_order_ignored");
        CHECK(result.diagnostics[0].line == 15);
        CHECK(result.diagnostics[1].code == "v3000_bond_order_ignored");
        CHECK(result.diagnostics[1].line == 17);
    }

    {
        std::ifstream input{fixture("corpus/sdf/heme/ideal.sdf")};
        auto reader = sdf::SdfReader{input, "HEM_ideal.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& hem_record = *result;
        CHECK(hem_record.molecule.atom_count() == 75);
        CHECK(hem_record.molecule.bond_count() == 80);
        REQUIRE(hem_record.diagnostics.size() == 1);
        CHECK(hem_record.diagnostics[0].code == "v3000_bond_order_ignored");
        CHECK(hem_record.diagnostics[0].line == 164);
        CHECK_FALSE(reader.next().has_value());
    }
}

TEST_CASE("native MOL2 input preserves names, charges, and bond orders", "[adapters][native]") {
    {
        std::ifstream input{fixture("synthetic/mol2/aromatic.mol2")};
        auto reader = mol2::Mol2Reader{input, "charged_aromatic.mol2"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& mol2_record = *result;
        CHECK(mol2_record.molecule.atom_count() == 3);
        CHECK(mol2_record.molecule.bond_count() == 2);
        CHECK(mol2_record.molecule.atom(0).name() == "N1");
        CHECK(mol2_record.molecule.atom(0).formal_charge() == 0);
        CHECK(mol2_record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(mol2_record.molecule.bond(1).order() == chargefw::core::BondOrder::DOUBLE);
        CHECK(mol2_record.diagnostics.size() == 1);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        auto input = std::istringstream{"@<TRIPOS>MOLECULE\namide\n3 2 0 0 0\nSMALL\nNO_CHARGES\n\n"
                                        "@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 O1 1.2 0 0 O.2\n"
                                        "3 N1 -1.3 0 0 N.am\n"
                                        "@<TRIPOS>BOND\n1 1 2 2\n2 1 3 am\n"};
        auto reader = mol2::Mol2Reader{input, "amide.mol2"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        CHECK(result->molecule.bond(0).order() == chargefw::core::BondOrder::DOUBLE);
        CHECK(result->molecule.bond(1).order() == chargefw::core::BondOrder::SINGLE);
    }
}

TEST_CASE("native SDF input preserves V2000, V3000, and multi-record streams",
          "[adapters][native]") {
    {
        std::ifstream input{fixture("synthetic/sdf/water.sdf")};
        auto reader = sdf::SdfReader{input, "water.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& water_record = *result;
        CHECK(water_record.molecule.atom_count() == 3);
        CHECK(water_record.molecule.bond_count() == 2);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/sdf/many_records.sdf")};
        auto reader = sdf::SdfReader{input, "set500.sdf"};
        std::size_t count = 0;

        while (const auto result = reader.next()) {
            ++count;
        }

        CHECK(count == 500);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/v3000.sdf")};
        auto reader = sdf::SdfReader{input, "set_v3000.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& v3000_record = *result;
        CHECK(v3000_record.molecule.atom_count() == 36);
        CHECK(v3000_record.molecule.bond_count() == 38);
        REQUIRE(v3000_record.import_metadata.has_value());
        CHECK(v3000_record.import_metadata->format == adapters::MolecularSourceFormat::sdf);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/mixed_v2000_v3000.sdf")};
        auto reader = sdf::SdfReader{input, "mixed_v2000_v3000.sdf"};
        const auto first = reader.next();
        const auto second = reader.next();
        REQUIRE(first.has_value());
        REQUIRE(second.has_value());
        CHECK(first->identity.record_id == "v2000");
        CHECK(second->identity.record_id == "v3000");
        CHECK(first->molecule.atom_count() == 1);
        CHECK(second->molecule.atom_count() == 1);
        CHECK_FALSE(reader.next().has_value());
    }
}

TEST_CASE("native readers reject malformed records", "[adapters][native]") {
    struct MalformedFixture {
        std::string_view path;
        void (*read_first)(std::istream& input);
    };
    for (const auto& [path, read_first] :
         {MalformedFixture{"synthetic/mol2/malformed_then_valid.mol2",
                           [](std::istream& input) {
                               static_cast<void>(mol2::Mol2Reader{input, "malformed.mol2"}.next());
                           }},
          MalformedFixture{"synthetic/mol2/missing_atom_section_then_valid.mol2",
                           [](std::istream& input) {
                               static_cast<void>(mol2::Mol2Reader{input, "malformed.mol2"}.next());
                           }},
          MalformedFixture{"synthetic/sdf/malformed_then_water.sdf", [](std::istream& input) {
                               static_cast<void>(sdf::SdfReader{input, "malformed.sdf"}.next());
                           }}}) {
        CAPTURE(path);
        std::ifstream input{fixture(path)};
        REQUIRE(input);
        CHECK_THROWS_AS(read_first(input), std::exception);
    }

    for (const auto malformed :
         {std::string_view{"bad-v2000\nchargefw\n\n  1  1  0  0  0  0  0  0  0  0999 V2000\n"
                           "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                           "  1  2  1  0  0  0  0\nM  END\n"},
          std::string_view{"bad-v3000\nchargefw\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                           "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 0 0 0 0\nM  V30 BEGIN ATOM\n"
                           "M  V30 1 C 0 0 0 0\nM  V30 1 O 1 0 0 0\nM  V30 END ATOM\n"
                           "M  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"}}) {
        std::istringstream input{std::string{malformed}};
        CHECK_THROWS_AS(mol::parse_mol(input, {}), std::exception);
    }

    for (const auto malformed : {
             std::string_view{
                 "bad-v2000-property\nchargefw\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                 "M  X\nM  END\n"},
             std::string_view{
                 "bad-v3000-bond\nchargefw\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                 "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n"
                 "M  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0\nM  V30 END ATOM\n"
                 "M  V30 BEGIN BOND\nM  V30 invalid 1 1 2\nM  V30 END BOND\n"
                 "M  V30 END CTAB\nM  END\n"},
             std::string_view{
                 "short-v3000-bond\nchargefw\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                 "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n"
                 "M  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0\nM  V30 END ATOM\n"
                 "M  V30 BEGIN BOND\nM  V30 1 1 1\nM  V30 END BOND\n"
                 "M  V30 END CTAB\nM  END\n"},
         }) {
        std::istringstream input{std::string{malformed}};
        CHECK_THROWS_AS(mol::parse_mol(input, {}), std::runtime_error);
    }
}
