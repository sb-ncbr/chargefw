#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/adapters/native/json_input.h>
#include <chargefw/adapters/native/mol2_output.h>
#include <chargefw/adapters/native/pqr_output.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <type_traits>

auto main() -> int {
    static_assert(std::is_class_v<chargefw::adapters::gemmi::mmcif_input::MmcifReader>);
    static_assert(std::is_class_v<chargefw::adapters::native::json_input::JsonReader>);
    static_assert(std::is_class_v<chargefw::adapters::gemmi::mmcif_output::MmcifWriter>);
    static_assert(std::is_class_v<chargefw::adapters::native::mol2_output::Mol2Writer>);
    static_assert(std::is_class_v<chargefw::adapters::native::pqr_output::PqrWriter>);

    return chargefw::parameters::load_default_parameter_sets().empty() ? 1 : 0;
}
