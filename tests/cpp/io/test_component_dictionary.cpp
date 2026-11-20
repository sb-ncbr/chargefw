#include "support/test_exceptions.h"

#include <chargefw/core/bond.h>
#include <chargefw/io/component_dictionary.h>
#include <chargefw/io/parse_error.h>
#include <snitch/snitch.hpp>

#include <filesystem>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

namespace io = chargefw::io;
using chargefw::core::BondOrder;

namespace {

// Formaldehyde, a zinc ion, and a hydroxide ion with a leaving hydrogen, in the format of
// chargefw-ccd. tests/fixtures/synthetic/components/components.txt.gz is this text compressed.
constexpr auto dictionary_text = "chargefw-ccd 1 2026.11.13.post1\n"
                                 "FOR 4 3\n"
                                 "C C 0\n"
                                 "O O 0\n"
                                 "H1 H 0\n"
                                 "H2 H 0\n"
                                 "0 1 2\n"
                                 "0 2 1\n"
                                 "0 3 1\n"
                                 "ION 1 0\n"
                                 "ZN Zn 0\n"
                                 "OH 2 1\n"
                                 "O O 0\n"
                                 "HO H 1\n"
                                 "0 1 1\n";

auto read(const std::string& text) -> io::ComponentDictionary {
    std::istringstream input{text};
    return io::ComponentDictionary{input, "components.txt"};
}

} // namespace

TEST_CASE("A component dictionary gives the atoms and bonds of a component by name",
          "[io][components]") {
    const auto dictionary = read(dictionary_text);

    CHECK(dictionary.version() == "2026.11.13.post1");
    CHECK(dictionary.size() == 3);

    // The last component, so that its place in the file matters.
    const auto hydroxide = dictionary.find("OH");
    REQUIRE(hydroxide.has_value());
    REQUIRE(hydroxide->atoms.size() == 2);
    CHECK(hydroxide->atoms[1].name == "HO");
    CHECK(hydroxide->atoms[1].atomic_number == 1);
    CHECK(hydroxide->atoms[1].leaving);
    CHECK_FALSE(hydroxide->atoms[0].leaving);

    const auto formaldehyde = dictionary.find("FOR");
    REQUIRE(formaldehyde.has_value());
    REQUIRE(formaldehyde->bonds.size() == 3);
    CHECK(formaldehyde->bonds[0].first == 0);
    CHECK(formaldehyde->bonds[0].second == 1);
    CHECK(formaldehyde->bonds[0].order == BondOrder::DOUBLE);
    CHECK(formaldehyde->bonds[2].second == 3);
    CHECK(formaldehyde->bonds[2].order == BondOrder::SINGLE);

    const auto ion = dictionary.find("ION");
    REQUIRE(ion.has_value());
    CHECK(ion->atoms[0].atomic_number == 30);
    CHECK(ion->bonds.empty());

    CHECK_FALSE(dictionary.find("ALA").has_value());
}

TEST_CASE("A component dictionary refuses a file it cannot read as a whole", "[io][components]") {
    struct Case {
        std::string text;
        io::ParseErrorKind kind;
        std::size_t line;
    };
    const auto cases = {
        // A later format.
        Case{"chargefw-ccd 2 2026.11.13\n", io::ParseErrorKind::unsupported, 1},
        // The file ends within a component.
        Case{"chargefw-ccd 1 2026.11.13\nOH 2 1\nO O 0\nHO H 1\n", io::ParseErrorKind::malformed,
             5},
        // A bond to an atom position that the component does not have.
        Case{"chargefw-ccd 1 2026.11.13\nOH 2 1\nO O 0\nHO H 1\n0 2 1\n",
             io::ParseErrorKind::malformed, 5},
        // Two atoms of one name, which a residue's atoms could not be matched to.
        Case{"chargefw-ccd 1 2026.11.13\nOH 2 1\nO O 0\nO H 1\n0 1 1\n",
             io::ParseErrorKind::malformed, 4},
    };
    for (const auto& [text, kind, line] : cases) {
        CAPTURE(text);
        const auto error = chargefw::test::thrown<io::ParseError>([&] { read(text); });
        REQUIRE(error.has_value());
        CHECK(error->kind() == kind);
        CHECK(error->location().source_name == "components.txt");
        CHECK(error->location().line == line);
    }

    // Another kind of file has no line to name: an mmCIF document, and a name that only starts
    // like the dictionary's.
    for (const auto* const text : {"data_1EHZ\n", "chargefw-ccd-other 1 2026.11.13\n"}) {
        CAPTURE(text);
        const auto error = chargefw::test::thrown<io::ParseError>([&] { read(text); });
        REQUIRE(error.has_value());
        CHECK(error->kind() == io::ParseErrorKind::no_records);
    }
}

TEST_CASE("A component dictionary is read from the compressed file of the package",
          "[io][components]") {
    const auto directory = std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" /
                           "synthetic" / "components";
    const auto dictionary = io::ComponentDictionary::read(directory / "components.txt.gz");

    CHECK(dictionary.version() == "2026.11.13.post1");
    const auto formaldehyde = dictionary.find("FOR");
    REQUIRE(formaldehyde.has_value());
    CHECK(formaldehyde->bonds.size() == 3);

    CHECK_THROWS_AS(io::ComponentDictionary::read(directory / "missing.txt.gz"),
                    std::runtime_error);

    // A file that does not decompress is refused by name, with no line.
    const auto corrupt = directory / "corrupt.txt.gz";
    const auto error = chargefw::test::thrown<io::ParseError>(
        [&] { static_cast<void>(io::ComponentDictionary::read(corrupt)); });
    REQUIRE(error.has_value());
    CHECK(error->kind() == io::ParseErrorKind::malformed);
    CHECK(error->location().source_name == corrupt.string());
    CHECK_FALSE(error->location().line.has_value());
}

TEST_CASE("A component dictionary is not taken from a stream that fails after a component",
          "[io][components]") {
    // Gives the text, which ends with a complete component, and then fails as a device would.
    class FailingBuffer : public std::streambuf {
      public:
        explicit FailingBuffer(std::string text) : text_{std::move(text)} {
            setg(text_.data(), text_.data(), text_.data() + text_.size());
        }

      protected:
        auto underflow() -> int_type override {
            throw std::runtime_error{"read error"};
        }

      private:
        std::string text_;
    };
    auto buffer = FailingBuffer{dictionary_text};
    auto input = std::istream{&buffer};

    // The text alone is a whole dictionary, so only the failure can stop it from being accepted.
    // The failure is not a refusal of the text: ParseError is a std::runtime_error too.
    auto failed_to_read = false;
    try {
        static_cast<void>(io::ComponentDictionary{input});
    } catch (const io::ParseError&) {
        failed_to_read = false;
    } catch (const std::runtime_error&) {
        failed_to_read = true;
    }
    CHECK(failed_to_read);
}
