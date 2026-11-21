#include <chargefw/io/component_dictionary.h>

#include <chargefw/core/periodic_table.h>
#include <chargefw/io/parse_error.h>

#include "parse_context.h"
#include "text_encoding.h"

#include <gemmi/gz.hpp>

#include <array>
#include <charconv>
#include <fstream>
#include <stdexcept>
#include <system_error>
#include <unordered_set>
#include <utility>

namespace chargefw::io {
namespace {

constexpr auto magic = std::string_view{"chargefw-ccd"};
constexpr auto format_version = std::string_view{"1"};

using Fields = std::array<std::string_view, 3>;

// Reads the dictionary text line by line and names the line of anything it refuses.
class Lines {
  public:
    Lines(const std::string_view text, const std::string& source_name, const std::size_t offset = 0)
        : text_{text}, source_name_{&source_name}, offset_{offset} {}

    [[nodiscard]] auto at_end() const noexcept -> bool {
        return offset_ >= text_.size();
    }
    [[nodiscard]] auto offset() const noexcept -> std::size_t {
        return offset_;
    }

    [[noreturn]] auto fail(const ParseErrorKind kind, std::string reason) const -> void {
        throw ParseError{kind, std::move(reason), {.source_name = *source_name_, .line = line_}};
    }

    // The three fields of the next line, which `what` describes.
    [[nodiscard]] auto next(const std::string_view what) -> Fields {
        ++line_;
        if (at_end()) {
            fail(ParseErrorKind::malformed,
                 "the file ends where " + std::string{what} + " belongs");
        }
        const auto end = std::min(text_.find('\n', offset_), text_.size());
        const auto line = text_.substr(offset_, end - offset_);
        offset_ = end + 1;

        const auto first = line.find(' ');
        const auto second = first == std::string_view::npos ? first : line.find(' ', first + 1);
        if (second == std::string_view::npos ||
            line.find(' ', second + 1) != std::string_view::npos) {
            fail(ParseErrorKind::malformed, std::string{what} + " does not have three fields");
        }
        const auto fields =
            Fields{line.substr(0, first), line.substr(first + 1, second - first - 1),
                   line.substr(second + 1)};
        for (const auto field : fields) {
            if (field.empty()) {
                fail(ParseErrorKind::malformed, std::string{what} + " has an empty field");
            }
        }
        return fields;
    }

    [[nodiscard]] auto number(const std::string_view field, const std::string_view what) const
        -> std::size_t {
        auto value = std::size_t{};
        const auto [end, error] = std::from_chars(field.data(), field.data() + field.size(), value);
        if (error != std::errc{} || end != field.data() + field.size()) {
            fail(ParseErrorKind::malformed,
                 std::string{what} + " '" + std::string{field} + "' is not a number");
        }
        return value;
    }

  private:
    std::string_view text_;
    const std::string* source_name_;
    std::size_t offset_;
    std::size_t line_ = 0;
};

[[nodiscard]] auto atomic_number(const Lines& lines, const std::string_view symbol) -> int {
    if (symbol == "X") {
        return 0; // The dictionary's own mark for an atom of no particular element.
    }
    if (symbol == "D") {
        return 1;
    }
    if (!core::periodic_table().contains(symbol)) {
        lines.fail(ParseErrorKind::malformed, "unknown element '" + std::string{symbol} + "'");
    }
    return core::periodic_table().element(symbol).atomic_number;
}

// The component whose first line has the given fields, with `lines` left after its last line.
[[nodiscard]] auto read_component(Lines& lines, const Fields& header) -> DictionaryComponent {
    const auto atom_count = lines.number(header[1], "atom count");
    const auto bond_count = lines.number(header[2], "bond count");

    auto component = DictionaryComponent{};
    auto names = std::unordered_set<std::string_view>{};
    for (std::size_t i = 0; i < atom_count; ++i) {
        const auto fields = lines.next("an atom");
        // Residue atoms are matched to these by name.
        if (!names.insert(fields[0]).second) {
            lines.fail(ParseErrorKind::malformed, "component " + std::string{header[0]} +
                                                      " lists atom " + std::string{fields[0]} +
                                                      " twice");
        }
        if (fields[2] != "0" && fields[2] != "1") {
            lines.fail(ParseErrorKind::malformed, "the leaving flag is not 0 or 1");
        }
        component.atoms.push_back({.name = std::string{fields[0]},
                                   .atomic_number = atomic_number(lines, fields[1]),
                                   .leaving = fields[2] == "1"});
    }
    for (std::size_t i = 0; i < bond_count; ++i) {
        const auto fields = lines.next("a bond");
        const auto first = lines.number(fields[0], "bond atom");
        const auto second = lines.number(fields[1], "bond atom");
        const auto order = lines.number(fields[2], "bond order");
        if (first >= atom_count || second >= atom_count) {
            lines.fail(ParseErrorKind::malformed,
                       "a bond of component " + std::string{header[0]} +
                           " names an atom position the component does not have");
        }
        if (first == second) {
            lines.fail(ParseErrorKind::malformed, "a bond joins an atom to itself");
        }
        if (order < 1 || order > 3) {
            lines.fail(ParseErrorKind::malformed, "the bond order is not 1, 2, or 3");
        }
        component.bonds.push_back({.first = first,
                                   .second = second,
                                   .order = core::bond_order_from_value(static_cast<int>(order))});
    }
    return component;
}

// The whole stream. Reading through the stream, not its buffer, records a failure in its state.
[[nodiscard]] auto contents(std::istream& input) -> std::string {
    auto text = std::string{};
    auto chunk = std::array<char, 65536>{};
    while (input.read(chunk.data(), static_cast<std::streamsize>(chunk.size())) ||
           input.gcount() > 0) {
        text.append(chunk.data(), static_cast<std::size_t>(input.gcount()));
    }
    detail::require_readable(input);
    return text;
}

} // namespace

auto ComponentDictionary::read(const std::filesystem::path& path) -> ComponentDictionary {
    const auto name = path.string();
    auto file = std::ifstream{path, std::ios::binary};
    if (!file) {
        throw std::runtime_error{"Unable to open component dictionary: " + name};
    }
    auto compressed = ::gemmi::MaybeGzipped{name};
    if (!compressed.is_compressed()) {
        return ComponentDictionary{file, name};
    }
    file.close();
    try {
        const auto buffer = detail::parsed([&] { return compressed.uncompress_into_buffer(); });
        return ComponentDictionary{std::string{buffer.data(), buffer.size()}, name};
    } catch (const ParseError&) {
        // A file that does not decompress is refused without a location; give it the file's name.
        detail::rethrow_located({.source_name = name});
    }
}

ComponentDictionary::ComponentDictionary(std::istream& input, const std::string& source_name)
    : ComponentDictionary{contents(input), source_name} {}

ComponentDictionary::ComponentDictionary(std::string text, const std::string& source_name)
    : text_{std::move(text)} {
    // The file is ASCII. Its version goes into results, which cannot hold other bytes as text.
    if (const auto line = detail::first_non_ascii_line(text_)) {
        throw ParseError{ParseErrorKind::invalid_encoding,
                         "a character that is not ASCII",
                         {.source_name = source_name, .line = line}};
    }
    auto lines = Lines{text_, source_name};
    // The first word decides whether this is a dictionary at all, before its line is judged.
    if (std::string_view{text_}.substr(0, text_.find_first_of(" \n")) != magic) {
        throw ParseError{ParseErrorKind::no_records,
                         "not a ChargeFW component dictionary",
                         {.source_name = source_name}};
    }
    const auto header = lines.next("the header");
    if (header[1] != format_version) {
        lines.fail(ParseErrorKind::unsupported,
                   "component dictionary format " + std::string{header[1]} +
                       "; this version of ChargeFW reads format " + std::string{format_version});
    }
    version_ = header[2];

    while (!lines.at_end()) {
        const auto offset = lines.offset();
        const auto component_header = lines.next("a component");
        if (!offsets_.emplace(component_header[0], offset).second) {
            lines.fail(ParseErrorKind::malformed,
                       "component " + std::string{component_header[0]} + " is listed twice");
        }
        // Checked here, so that find() cannot fail later.
        static_cast<void>(read_component(lines, component_header));
    }
}

auto ComponentDictionary::version() const noexcept -> const std::string& {
    return version_;
}

auto ComponentDictionary::size() const noexcept -> std::size_t {
    return offsets_.size();
}

auto ComponentDictionary::find(const std::string_view name) const
    -> std::optional<DictionaryComponent> {
    const auto found = offsets_.find(std::string{name});
    if (found == offsets_.end()) {
        return std::nullopt;
    }
    static const auto no_source = std::string{};
    auto lines = Lines{text_, no_source, found->second};
    const auto header = lines.next("a component");
    return read_component(lines, header);
}

} // namespace chargefw::io
