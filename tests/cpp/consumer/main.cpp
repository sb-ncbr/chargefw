// The umbrella header, so that the names of all readers and writers are checked together.
#include <chargefw/io/all.h>
#include <chargefw/parameters/io/parameter_set_io.h>
#include <chargefw/version.h>

#include <fstream>
#include <sstream>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

template <typename Reader> auto check_reader_moves(Reader reader) -> bool {
    auto moved = std::move(reader);
    const auto record = moved.next();
    reader = std::move(moved);
    return record && record->molecule.atom_count() > 0 &&
           record->identity.source_name == "consumer" && record->identity.record_index == 0 &&
           !reader.next();
}

} // namespace

static_assert(CHARGEFW_VERSION_MAJOR == EXPECTED_VERSION_MAJOR &&
              CHARGEFW_VERSION_MINOR == EXPECTED_VERSION_MINOR &&
              CHARGEFW_VERSION_PATCH == EXPECTED_VERSION_PATCH);
static_assert(std::string_view{CHARGEFW_VERSION_STRING}.starts_with(EXPECTED_VERSION));

auto main(int argc, char* argv[]) -> int {
    static_assert(std::is_class_v<chargefw::io::MmcifWriter>);
    static_assert(std::is_class_v<chargefw::io::Mol2Writer>);
    static_assert(std::is_class_v<chargefw::io::PqrWriter>);

    if (argc != 2) {
        return 1;
    }
    auto json_input = std::istringstream{R"({"schema_version":"1.0","molecules":[
        {"atoms":[{"atomic_number":1,"formal_charge":0}]}]})"};
    auto mmcif_input = std::ifstream{argv[1]};
    if (!check_reader_moves(chargefw::io::JsonReader{json_input, "consumer"}) ||
        !check_reader_moves(chargefw::io::MmcifReader{mmcif_input, "consumer"})) {
        return 1;
    }

    return chargefw::parameters::load_default_parameter_sets().empty() ? 1 : 0;
}
