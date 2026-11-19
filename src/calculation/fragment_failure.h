#pragma once

#include <chargefw/methods/method_failure.h>

#include <cstddef>

namespace chargefw::calculation::detail {

// A method failure in the fragment around an atom, numbered as in the input molecule.
class FragmentFailure : public methods::MethodFailure {
  public:
    FragmentFailure(const methods::MethodFailure& failure, const std::size_t atom_index)
        : methods::MethodFailure{failure}, atom_index_{atom_index} {}

    [[nodiscard]] auto atom_index() const noexcept -> std::size_t {
        return atom_index_;
    }

  private:
    std::size_t atom_index_;
};

} // namespace chargefw::calculation::detail
