#pragma once

#include <chargefw/core/bond.h>

#include <cstddef>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace chargefw::io {

struct DictionaryAtom {
    std::string name;
    int atomic_number = 0; // Zero when the dictionary gives no element.
    // Whether the atom is absent when the component is linked in a polymer, such as OXT.
    bool leaving = false;
};

// A bond between two atoms of a component, given by their positions in its atom list.
struct DictionaryBond {
    std::size_t first = 0;
    std::size_t second = 0;
    core::BondOrder order = core::BondOrder::SINGLE;
};

struct DictionaryComponent {
    std::vector<DictionaryAtom> atoms;
    std::vector<DictionaryBond> bonds;
};

// The atoms and bonds of chemical components by residue name, read from the data file of the
// chargefw-ccd package. Reading checks the whole file and throws ParseError when it is not a
// dictionary this version reads.
class ComponentDictionary {
  public:
    // Reads the data file, which is gzip-compressed when its name ends in .gz, as the package's is.
    // Throws std::runtime_error for a file that cannot be opened.
    [[nodiscard]] static auto read(const std::filesystem::path& path) -> ComponentDictionary;

    // Reads the decompressed text of the data file.
    explicit ComponentDictionary(std::istream& input, const std::string& source_name = {});

    // The version of the package the file came from, such as "2026.11.13" or "2026.11.13.post1".
    [[nodiscard]] auto version() const noexcept -> const std::string&;

    // The number of components.
    [[nodiscard]] auto size() const noexcept -> std::size_t;

    // The component with the given name, or nothing when the dictionary has none.
    [[nodiscard]] auto find(std::string_view name) const -> std::optional<DictionaryComponent>;

  private:
    ComponentDictionary(std::string text, const std::string& source_name);

    std::string text_;
    std::string version_;
    // Where each component starts in text_. Components are parsed when asked for.
    std::unordered_map<std::string, std::size_t> offsets_;
};

} // namespace chargefw::io
