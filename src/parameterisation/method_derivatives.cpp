#include "parameterisation/method_derivatives.h"

#include "methods/builtin/eem.h"
#include "methods/builtin/sqe.h"
#include "methods/builtin/sqeq0.h"
#include "methods/builtin/sqeqp.h"

#include <chargefw/methods/method.h>

namespace chargefw::parameterisation::detail {

auto derivatives_of(const methods::Method& method) -> const MethodDerivatives* {
    // The built-in method classes are final, so each match is exactly that method.
    if (dynamic_cast<const methods::builtin::EEMMethod*>(&method) != nullptr) {
        return &eem_derivatives();
    }
    if (dynamic_cast<const methods::builtin::SQEMethod*>(&method) != nullptr) {
        return &sqe_derivatives();
    }
    if (dynamic_cast<const methods::builtin::SQEq0Method*>(&method) != nullptr) {
        return &sqeq0_derivatives();
    }
    if (dynamic_cast<const methods::builtin::SQEqpMethod*>(&method) != nullptr) {
        return &sqeqp_derivatives();
    }
    return nullptr;
}

} // namespace chargefw::parameterisation::detail
