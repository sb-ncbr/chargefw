#include <chargefw/io/conformer_selection.h>

#include <stdexcept>
#include <string>

namespace chargefw::io {

auto conformer_selection_from_string(const std::string_view value) -> ConformerSelection {
    if (value == "first") {
        return ConformerSelection::first;
    }
    if (value == "all") {
        return ConformerSelection::all;
    }
    throw std::invalid_argument{"unknown conformer selection: " + std::string{value}};
}

auto to_string(const ConformerSelection value) -> std::string_view {
    switch (value) {
    case ConformerSelection::first:
        return "first";
    case ConformerSelection::all:
        return "all";
    }
    throw std::invalid_argument{"unknown conformer selection"};
}

} // namespace chargefw::io
