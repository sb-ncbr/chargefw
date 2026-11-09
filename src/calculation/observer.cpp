#include <chargefw/calculation/observer.h>

#include <stdexcept>
#include <string_view>

namespace chargefw::calculation {

auto to_string(const CalculationPhase value) -> std::string_view {
    switch (value) {
    case CalculationPhase::computation_started:
        return "computation_started";
    case CalculationPhase::computation_finished:
        return "computation_finished";
    case CalculationPhase::target_started:
        return "target_started";
    case CalculationPhase::target_finished:
        return "target_finished";
    case CalculationPhase::fragment_progress:
        return "fragment_progress";
    }
    throw std::invalid_argument{"unknown calculation phase"};
}

} // namespace chargefw::calculation
