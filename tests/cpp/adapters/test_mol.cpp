#include <algorithm>
#include <chargefw/adapters/native/mol2_input.h>
#include <chargefw/adapters/native/mol_input.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/core/bond.h>
#include <exception>
#include <filesystem>
#include <fstream>
#include <istream>
#include <iterator>
#include <locale>
#include <optional>
#include <snitch/snitch.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

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
        // Atom-block charge codes count unless M  CHG or M  RAD lines replace them all.
        constexpr std::string_view atoms =
            "\nchargefw\n\n  3  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n"
            "    1.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n"
            "    2.0000    0.0000    0.0000 Fe  0  2\n";
        std::istringstream block_only{std::string{atoms} + "M  END\n"};
        const auto from_block = mol::parse_mol(block_only, {});
        CHECK(from_block.molecule.atom(0).formal_charge() == 1);
        CHECK(from_block.molecule.atom(1).formal_charge() == -1);
        CHECK(from_block.molecule.atom(2).formal_charge() == 2);

        std::istringstream with_chg{std::string{atoms} + "M  CHG  1   3   3\nM  END\n"};
        const auto from_chg = mol::parse_mol(with_chg, {});
        CHECK(from_chg.molecule.atom(0).formal_charge() == 0);
        CHECK(from_chg.molecule.atom(1).formal_charge() == 0);
        CHECK(from_chg.molecule.atom(2).formal_charge() == 3);

        std::istringstream with_rad{std::string{atoms} + "M  RAD  1   3   2\nM  END\n"};
        const auto from_rad = mol::parse_mol(with_rad, {});
        for (std::size_t index = 0; index < 3; ++index) {
            CHECK(from_rad.molecule.atom(index).formal_charge() == 0);
        }
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
        CHECK(result.issues.empty());
    }
}

TEST_CASE("native MOL and SDF input ignore unused properties and report left-out bond orders",
          "[adapters][native]") {
    {
        std::ifstream input{fixture("synthetic/mol/v2000/ignored_properties.mol")};
        const auto result = mol::parse_mol(input, {});
        CHECK(result.molecule.atom_count() == 2);
        CHECK(result.molecule.bond_count() == 1);
        CHECK(result.molecule.atom(0).formal_charge() == 1);
        CHECK(result.molecule.atom(1).formal_charge() == -1);
        CHECK(result.issues.empty());
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
        REQUIRE(result.issues.size() == 2);
        CHECK(result.issues[0].kind ==
              chargefw::adapters::ImportIssueKind::v3000_bond_order_ignored);
        CHECK(result.issues[0].line == 15);
        CHECK(result.issues[1].kind ==
              chargefw::adapters::ImportIssueKind::v3000_bond_order_ignored);
        CHECK(result.issues[1].line == 17);
    }

    {
        std::ifstream input{fixture("corpus/sdf/heme/ideal.sdf")};
        auto reader = sdf::SdfReader{input, "HEM_ideal.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& hem_record = *result;
        CHECK(hem_record.molecule.atom_count() == 75);
        CHECK(hem_record.molecule.bond_count() == 80);
        REQUIRE(hem_record.issues.size() == 1);
        CHECK(hem_record.issues[0].kind ==
              chargefw::adapters::ImportIssueKind::v3000_bond_order_ignored);
        CHECK(hem_record.issues[0].line == 164);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::istringstream aromatic{
            "aromatic\nchargefw\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0\n"
            "    1.4000    0.0000    0.0000 C   0  0\n  1  2  4\nM  END\n"};
        const auto result = mol::parse_mol(aromatic, {});
        CHECK(result.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        REQUIRE(result.issues.size() == 1);
        CHECK(result.issues[0].kind ==
              chargefw::adapters::ImportIssueKind::aromatic_bonds_converted);
    }
}

TEST_CASE("native reader diagnostics use file lines across records", "[adapters][native]") {
    std::istringstream unnamed_input{"broken\n"};
    CHECK_THROWS_MATCHES(mol::parse_mol(unnamed_input, {}), std::runtime_error,
                         snitch::matchers::with_what_contains{"line 2: record 1 ('broken'):"});

    for (const auto newline : {std::string_view{"\n"}, std::string_view{"\r\n"}}) {
        const auto input_text = [newline](const std::string_view text) {
            std::string result;
            for (const char character : text) {
                if (character == '\n') {
                    result += newline;
                } else {
                    result += character;
                }
            }
            return result;
        };
        const std::string mol_record{
            "\n                    2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0\nM  ISO  1   1  13\nM  END\n"
            "> <comment>\nignored\n\n$$$$\n"};
        std::istringstream sdf_input{input_text(mol_record + mol_record + "broken\n")};
        auto sdf_reader = sdf::SdfReader{sdf_input, "input.sdf"};
        REQUIRE(sdf_reader.next().has_value());
        const auto second = sdf_reader.next();
        REQUIRE(second.has_value());
        REQUIRE(second->issues.size() == 1);
        CHECK(second->issues[0].line == 13);
        CHECK_THROWS_MATCHES(sdf_reader.next(), std::runtime_error,
                             snitch::matchers::with_what_contains{"input.sdf:24: record 3"});

        const std::string mol2_record{"@<TRIPOS>MOLECULE\nligand\n1 0\n@<TRIPOS>ATOM\n"
                                      "1 C1 0 0 0 C.3 1 LIG 0.2\n@<TRIPOS>BOND\n"};
        std::istringstream mol2_input{
            input_text("# preamble\n" + mol2_record + mol2_record +
                       "@<TRIPOS>MOLECULE\nbroken\n1 1\n@<TRIPOS>ATOM\n"
                       "1 C1 0 0 0 C.3\n@<TRIPOS>BOND\n1 1 9 1\n@<TRIPOS>SUBSTRUCTURE\nignored\n")};
        auto mol2_reader = mol2::Mol2Reader{mol2_input, "input.mol2"};
        REQUIRE(mol2_reader.next().has_value());
        const auto second_mol2 = mol2_reader.next();
        REQUIRE(second_mol2.has_value());
        REQUIRE(second_mol2->issues.size() == 1);
        CHECK(second_mol2->issues[0].line == 12);
        CHECK_THROWS_MATCHES(mol2_reader.next(), std::runtime_error,
                             snitch::matchers::with_what_contains{"input.mol2:20: record 3"});
    }
}

TEST_CASE("native readers reject non-ASCII text on any line", "[adapters][native]") {
    const std::string mol_body{"\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                               "    0.0000    0.0000    0.0000 C   0  0\nM  END\n"};
    // A Latin-1 title, and UTF-8 text in a data item, which ChargeFW does not otherwise read.
    std::istringstream title_input{"Caf\xe9" + mol_body + "$$$$\n"};
    auto title_reader = sdf::SdfReader{title_input, "input.sdf"};
    CHECK_THROWS_MATCHES(title_reader.next(), std::runtime_error,
                         snitch::matchers::with_what_contains{"input.sdf:1: record 1:"});
    std::istringstream data_input{"methane" + mol_body + "> <name>\n\xce\xb1-pinene\n\n$$$$\n"};
    auto data_reader = sdf::SdfReader{data_input, "input.sdf"};
    CHECK_THROWS_MATCHES(
        data_reader.next(), std::runtime_error,
        snitch::matchers::with_what_contains{"input.sdf:8: record 1 ('methane'):"});
    std::istringstream mol_input{"methane" + mol_body + "Caf\xe9\n"};
    auto mol_reader = mol::MolReader{mol_input, "input.mol"};
    CHECK_THROWS_MATCHES(
        mol_reader.next(), std::runtime_error,
        snitch::matchers::with_what_contains{"input.mol:7: record 1 ('methane'):"});

    std::istringstream mol2_input{"# \xe5\x88\x86\xe5\xad\x90\n@<TRIPOS>MOLECULE\nligand\n1 0\n"
                                  "@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n@<TRIPOS>BOND\n"};
    auto mol2_reader = mol2::Mol2Reader{mol2_input, "input.mol2"};
    CHECK_THROWS_MATCHES(mol2_reader.next(), std::runtime_error,
                         snitch::matchers::with_what_contains{"input.mol2:1: record 1:"});
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
        REQUIRE(mol2_record.issues.size() == 2);
        CHECK(mol2_record.issues[0].kind ==
              chargefw::adapters::ImportIssueKind::partial_charges_ignored);
        CHECK(mol2_record.issues[1].kind ==
              chargefw::adapters::ImportIssueKind::aromatic_bonds_converted);
        CHECK_FALSE(reader.next().has_value());
    }

    // Aromatic bond types are read as single bonds with a warning; an amide C-N bond is single.
    for (const auto& [type, aromatic] :
         std::vector<std::pair<std::string, bool>>{{"am", false}, {"ar", true}, {"4", true}}) {
        auto input = std::istringstream{"@<TRIPOS>MOLECULE\nbonds\n3 2 0 0 0\nSMALL\nNO_CHARGES\n\n"
                                        "@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 O1 1.2 0 0 O.2\n"
                                        "3 N1 -1.3 0 0 N.am\n"
                                        "@<TRIPOS>BOND\n1 1 2 2\n2 1 3 " +
                                        type + "\n"};
        auto reader = mol2::Mol2Reader{input, "bonds.mol2"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        CHECK(result->molecule.bond(0).order() == chargefw::core::BondOrder::DOUBLE);
        CHECK(result->molecule.bond(1).order() == chargefw::core::BondOrder::SINGLE);
        REQUIRE(result->issues.size() == (aromatic ? 1U : 0U));
        if (aromatic) {
            CHECK(result->issues[0].kind ==
                  chargefw::adapters::ImportIssueKind::aromatic_bonds_converted);
            CHECK(result->issues[0].location.bond_indices == std::vector<std::size_t>{1});
        }
    }
}

TEST_CASE("MOL2 markers tolerate trailing whitespace without losing records",
          "[adapters][native]") {
    std::istringstream input{"@<TRIPOS>MOLECULE \t\r\nfirst\n1 0\nSMALL\nNO_CHARGES\n"
                             "@<TRIPOS>ATOM \t\n1 Na 0 0 0 Na\n \t\n# atoms end\n"
                             "@<TRIPOS>BOND \t\n \t# no bonds\n"
                             "@<TRIPOS>UNITY_ATOM_ATTR \t\n1 1\ncharge 1\n"
                             "@<TRIPOS>MOLECULE \t\nsecond\n1 0\nSMALL\nNO_CHARGES\n"
                             "@<TRIPOS>ATOM \t\n1 Cl 1 0 0 Cl\n"
                             "@<TRIPOS>UNITY_ATOM_ATTR \t\n1 1\ncharge -1\n"
                             "@<TRIPOS>BOND \t\n"};
    mol2::Mol2Reader reader{input, "padded.mol2"};
    const auto first = reader.next();
    REQUIRE(first.has_value());
    CHECK(first->molecule.name() == "first");
    CHECK(first->molecule.atom(0).atomic_number() == 11);
    CHECK(first->molecule.atom(0).formal_charge() == 1);
    const auto second = reader.next();
    REQUIRE(second.has_value());
    CHECK(second->molecule.name() == "second");
    CHECK(second->molecule.atom(0).atomic_number() == 17);
    CHECK(second->molecule.atom(0).formal_charge() == -1);
    CHECK_FALSE(reader.next().has_value());
}

TEST_CASE("MOL2 rejects repeated atom and bond sections", "[adapters][native]") {
    for (const auto repeated : {"@<TRIPOS>ATOM \t\n1 Cl 1 0 0 Cl\n", "@<TRIPOS>BOND \t\n"}) {
        CAPTURE(repeated);
        std::istringstream input{std::string{"@<TRIPOS>MOLECULE\nion\n1 0\nSMALL\nNO_CHARGES\n"
                                             "@<TRIPOS>ATOM\n1 Na 0 0 0 Na\n@<TRIPOS>BOND\n"} +
                                 repeated};
        mol2::Mol2Reader reader{input, "duplicate.mol2"};
        CHECK_THROWS_AS(reader.next(), std::runtime_error);
    }
}

TEST_CASE("MOL2 rejects atom and bond counts smaller than their sections", "[adapters][native]") {
    for (const auto counts : {"1 0", "2 0", "2 1"}) {
        CAPTURE(counts);
        std::istringstream input{
            std::string{"@<TRIPOS>MOLECULE\npair\n"} + counts +
            "\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n2 H 1 0 0 H\n@<TRIPOS>BOND\n1 1 2 1\n"};
        mol2::Mol2Reader reader{input, "counts.mol2"};
        if (std::string_view{counts} == "2 1") {
            const auto record = reader.next();
            REQUIRE(record.has_value());
            CHECK(record->molecule.atom_count() == 2);
            CHECK(record->molecule.bond_count() == 1);
            CHECK_FALSE(reader.next().has_value());
        } else {
            CHECK_THROWS_AS(reader.next(), std::runtime_error);
        }
    }
}

TEST_CASE("native MOL2 input reads formal charges from UNITY_ATOM_ATTR", "[adapters][native]") {
    // Open Babel writes UNITY_ATOM_ATTR between ATOM and BOND; the format allows any section order,
    // and the second record has it after BOND, next to an attribute that is not a charge.
    std::istringstream input{"@<TRIPOS>MOLECULE\nformate\n3 2 0 0 0\nSMALL\nNO_CHARGES\n\n"
                             "@<TRIPOS>ATOM\n7 C1 0 0 0 C.2\n9 O1 1.2 0 0 O.co2\n"
                             "11 O2 -1.2 0 0 O.co2\n"
                             "@<TRIPOS>UNITY_ATOM_ATTR\n11 1\ncharge -1\n"
                             "@<TRIPOS>BOND\n1 7 9 2\n2 7 11 1\n"
                             "@<TRIPOS>SUBSTRUCTURE\n1 FOR 1\n"
                             "@<TRIPOS>MOLECULE\nammonium\n1 0 0 0 0\nSMALL\nNO_CHARGES\n\n"
                             "@<TRIPOS>ATOM\n1 N1 0 0 0 N.4\n@<TRIPOS>BOND\n"
                             "@<TRIPOS>UNITY_ATOM_ATTR\n1 2\nS N\ncharge 1\n"};
    auto reader = mol2::Mol2Reader{input, "charged.mol2"};

    const auto formate = reader.next();
    REQUIRE(formate.has_value());
    CHECK(formate->molecule.atom(0).formal_charge() == 0);
    CHECK(formate->molecule.atom(1).formal_charge() == 0);
    CHECK(formate->molecule.atom(2).formal_charge() == -1);
    CHECK(formate->molecule.bond_count() == 2);
    CHECK(formate->issues.empty());

    const auto ammonium = reader.next();
    REQUIRE(ammonium.has_value());
    CHECK(ammonium->identity.record_id == "ammonium");
    CHECK(ammonium->molecule.atom(0).formal_charge() == 1);
    CHECK_FALSE(reader.next().has_value());

    // An attribute count must match the attribute lines. A count of 6 would otherwise take the next
    // record's header as attributes and its atoms as those of this record.
    constexpr std::string_view ammonium_record =
        "@<TRIPOS>MOLECULE\nammonium\n1 0 0 0 0\nSMALL\nNO_CHARGES\n\n"
        "@<TRIPOS>ATOM\n1 N1 0 0 0 N.4\n@<TRIPOS>BOND\n";
    for (const auto attributes : {"1 6\ncharge 1\n", "1 -1\n", "1 2\ncharge 1\n\n"}) {
        CAPTURE(attributes);
        std::istringstream malformed{std::string{ammonium_record} + "@<TRIPOS>UNITY_ATOM_ATTR\n" +
                                     attributes + std::string{ammonium_record}};
        auto malformed_reader = mol2::Mol2Reader{malformed, "malformed.mol2"};
        CHECK_THROWS_AS(static_cast<void>(malformed_reader.next()), std::runtime_error);
    }
}

TEST_CASE("native SDF input preserves V2000, V3000, and multi-record streams",
          "[adapters][native]") {
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

TEST_CASE("native SDF input imports declared 2D drawings without coordinates",
          "[adapters][native]") {
    const auto read = [](const std::string_view name) {
        std::ifstream input{fixture("synthetic/sdf/" + std::string{name})};
        auto reader = sdf::SdfReader{input, std::string{name}};
        auto record = reader.next();
        REQUIRE(record.has_value());
        CHECK_FALSE(reader.next().has_value());
        return std::move(*record);
    };

    // The same water, once as a PubChem 2D drawing and once with 3D coordinates.
    const auto drawing = read("water.sdf");
    CHECK(drawing.molecule.atom_count() == 3);
    CHECK(drawing.molecule.bond_count() == 2);
    CHECK(drawing.molecule.conformer_count() == 0);
    REQUIRE(drawing.import_metadata.has_value());
    CHECK(drawing.import_metadata->conformers.empty());
    REQUIRE(drawing.issues.size() == 1);
    CHECK(drawing.issues.front().kind ==
          chargefw::adapters::ImportIssueKind::coordinates_2d_ignored);

    const auto geometry = read("water_3d.sdf");
    CHECK(geometry.molecule.conformer_count() == 1);
    CHECK(geometry.issues.empty());

    // Only the header declaration counts: V3000 records too, and flat records without it are kept.
    auto v3000_stream = std::ifstream{fixture("synthetic/mol/v3000/charged_atoms.mol")};
    const auto v3000 = std::string{std::istreambuf_iterator<char>{v3000_stream}, {}};
    auto v3000_drawing = v3000;
    v3000_drawing.replace(v3000_drawing.find("  ChargeFW\n"), 11, "  ChargeFW          2D\n");
    struct Case {
        std::string_view description;
        std::string text;
        std::size_t conformers;
    };
    for (const auto& [description, text, conformers] : std::vector<Case>{
             {"V3000 without a declaration", v3000, 1},
             {"V3000 declared 2D", v3000_drawing, 0},
             {"flat V2000 without a declaration",
              "water\n  ChargeFW\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
              "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
              "    0.9570    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
              "   -0.2400    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
              "  1  2  1  0\n  1  3  1  0\nM  END\n",
              1}}) {
        CAPTURE(description);
        auto input = std::istringstream{text};
        const auto record = mol::parse_mol(input, {.source = "input.mol"});
        CHECK(record.molecule.conformer_count() == conformers);
        CHECK(std::ranges::count(
                  record.issues, chargefw::adapters::ImportIssueKind::coordinates_2d_ignored,
                  &chargefw::adapters::ImportIssue::kind) == (conformers == 0 ? 1 : 0));
    }
}

TEST_CASE("V3000 zero-bond records accept omitted and empty bond blocks", "[adapters][native]") {
    for (const auto bond_block : {"", "M  V30 BEGIN BOND\nM  V30 END BOND\n"}) {
        CAPTURE(bond_block);
        // RDKit omits the block; Open Babel writes an empty block for the same sodium ion.
        std::istringstream input{
            std::string{"sodium\nchargefw\n\n"
                        "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                        "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\n"
                        "M  V30 1 Na 0 0 0 0 CHG=1\nM  V30 END ATOM\n"} +
            bond_block + "M  V30 END CTAB\nM  END\n"};
        const auto record = mol::parse_mol(input, {});
        REQUIRE(record.molecule.atom_count() == 1);
        CHECK(record.molecule.atom(0).atomic_number() == 11);
        CHECK(record.molecule.atom(0).formal_charge() == 1);
        CHECK(record.molecule.bond_count() == 0);
        REQUIRE(record.import_metadata.has_value());
        CHECK(record.import_metadata->source_connectivity ==
              adapters::SourceConnectivity::explicitly_empty);
    }
}

TEST_CASE("V3000 stereo collections preserve the graph and SDF record boundaries",
          "[adapters][native]") {
    // RDKit output for C[C@H](F)[C@H](Cl)Br |&1:1,3|, with a continued collection record.
    const auto stereo =
        std::string{"stereo\n     RDKit          2D\n\n"
                    "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                    "M  V30 BEGIN CTAB\nM  V30 COUNTS 6 5 0 0 0\nM  V30 BEGIN ATOM\n"
                    "M  V30 1 C 0.000000 0.000000 0.000000 0\n"
                    "M  V30 2 C 1.299038 0.750000 0.000000 0\n"
                    "M  V30 3 F 1.299038 2.250000 0.000000 0\n"
                    "M  V30 4 C 2.598076 -0.000000 0.000000 0\n"
                    "M  V30 5 Cl 3.897114 0.750000 0.000000 0\n"
                    "M  V30 6 Br 2.598076 -1.500000 0.000000 0\n"
                    "M  V30 END ATOM\nM  V30 BEGIN BOND\n"
                    "M  V30 1 1 2 1 CFG=3\nM  V30 2 1 2 3\nM  V30 3 1 2 4\n"
                    "M  V30 4 1 4 5 CFG=3\nM  V30 5 1 4 6\nM  V30 END BOND\n"
                    "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STERAC1 ATOMS=(2 2 -\n"
                    "M  V30 4)\nM  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n"};
    std::istringstream input{stereo + "$$$$\n" + v3000_with_coordinates("1 2 3") + "$$$$\n"};
    sdf::SdfReader reader{input, "stereo.sdf"};
    const auto first = reader.next();
    REQUIRE(first.has_value());
    REQUIRE(first->molecule.atom_count() == 6);
    REQUIRE(first->molecule.bond_count() == 5);
    CHECK(first->molecule.atom(5).atomic_number() == 35);
    CHECK(first->molecule.bond(0).first_atom_index() == 1);
    CHECK(first->molecule.bond(0).second_atom_index() == 0);
    const auto second = reader.next();
    REQUIRE(second.has_value());
    CHECK(second->molecule.atom_count() == 1);
    CHECK(second->molecule.conformer(0).positions()[0].z == 3.0);
    CHECK_FALSE(reader.next().has_value());
}

TEST_CASE("MOL input reads atoms, bonds, and charges and ignores the rest", "[adapters][native]") {
    // An atom alias and a group abbreviation, whose texts look like property lines, a polymer
    // S-group, and an isotope around the charge line.
    std::istringstream v2000{"ignored\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
                             "    0.0000    0.0000    0.0000 C   0  0\n"
                             "    1.4000    0.0000    0.0000 O   0  0\n  1  2  1  0\n"
                             "A    1\nM  END\nG    1  2\nM  CHG  1   1  -1\n"
                             "M  STY  1   1 SRU\nM  SAL   1  2   1   2\n"
                             "M  CHG  1   2  -1\nM  ISO  1   1  13\nM  END\n"};
    const auto v2000_record = mol::parse_mol(v2000, {});
    CHECK(v2000_record.molecule.atom_count() == 2);
    CHECK(v2000_record.molecule.bond_count() == 1);
    CHECK(v2000_record.molecule.atom(0).formal_charge() == 0);
    CHECK(v2000_record.molecule.atom(1).formal_charge() == -1);
    CHECK(v2000_record.issues.empty());

    const auto source = v3000_with_coordinates("0 0 0");
    const auto atoms_and_bonds = source.substr(0, source.find("M  V30 END CTAB\n"));
    for (const auto tail :
         {"M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(1 1) LABEL=Me\nM  V30 END SGROUP\n"
          "M  V30 END CTAB\nM  END\n",
          "M  V30 BEGIN OBJ3D\nM  V30 1 POINT\nM  V30 END OBJ3D\nM  V30 END CTAB\nM  END\n",
          "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION\n"
          "M  V30 END CTAB\nM  END\n",
          "M  V30 END CTAB\nM  V30 BEGIN RGROUP 1\nM  V30 END RGROUP\nM  END\n"}) {
        CAPTURE(tail);
        std::istringstream input{atoms_and_bonds + tail};
        const auto record = mol::parse_mol(input, {});
        CHECK(record.molecule.atom_count() == 1);
        CHECK(record.issues.empty());
    }

    // Ignoring stops at M  END; without it, the SDF record delimiter or the end of input is an
    // error.
    const std::string v2000_atom{"\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                                 "    0.0000    0.0000    0.0000 C   0  0\n"};
    for (const auto& text :
         {"no-end" + v2000_atom + "> <name>\nvalue\n\n$$$$\n",
          // A record that ends in place of alias text must not take the next record's charge.
          "alias-end" + v2000_atom + "A    1\n$$$$\nnext" + v2000_atom +
              "M  CHG  1   1  -1\nM  END\n$$$$\n",
          atoms_and_bonds + "M  V30 END CTAB\n$$$$\n", atoms_and_bonds + "M  V30 BEGIN SGROUP\n"}) {
        std::istringstream input{text};
        CHECK_THROWS_AS(mol::parse_mol(input, {}), std::runtime_error);
    }

    // The atoms and bonds that are read must match their counts and come in one block each.
    const auto atoms = source.substr(0, source.find("M  V30 BEGIN BOND\n"));
    auto missing_bond = atoms + "M  V30 END CTAB\nM  END\n";
    missing_bond.replace(missing_bond.find("COUNTS 1 0"), std::string_view{"COUNTS 1 0"}.size(),
                         "COUNTS 1 1");
    for (const auto& text :
         {missing_bond,
          atoms + "M  V30 BEGIN BOND\nM  V30 1 1 1 1\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n",
          atoms_and_bonds + "M  V30 BEGIN BOND\nM  V30 1 1 1 1\nM  V30 END BOND\n"
                            "M  V30 END CTAB\nM  END\n",
          atoms_and_bonds + "M  V30 BEGIN ATOM\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\n"
                            "M  V30 END CTAB\nM  END\n"}) {
        std::istringstream input{text};
        CHECK_THROWS_AS(mol::parse_mol(input, {}), std::runtime_error);
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
                 "bad-v2000-charge\nchargefw\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.0000    0.0000    0.0000 C   0  8  0  0  0  0  0  0  0  0  0  0\n"
                 "M  END\n"},
             std::string_view{
                 "unknown-element\nchargefw\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.0000    0.0000    0.0000 Xx  0  0  0  0  0  0  0  0  0  0  0  0\n"
                 "M  END\n"},
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

    std::istringstream unknown_mol2_element{
        "@<TRIPOS>MOLECULE\nunknown\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 X 0 0 0 Xx\n"};
    mol2::Mol2Reader mol2_reader{unknown_mol2_element, "unknown.mol2"};
    CHECK_THROWS_AS(mol2_reader.next(), std::runtime_error);
}
