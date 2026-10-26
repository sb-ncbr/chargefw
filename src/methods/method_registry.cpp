#include <chargefw/methods/method_registry.h>

#include "methods/builtin_methods.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::methods {
namespace {

auto validate_methods(const std::vector<std::unique_ptr<Method>>& methods) -> void {
    for (const auto& method : methods) {
        if (method == nullptr) {
            throw std::invalid_argument{"method registry contains null method"};
        }

        if (method->id().empty()) {
            throw std::invalid_argument{"method registry contains method with empty id"};
        }
    }

    for (auto first = methods.begin(); first != methods.end(); ++first) {
        const auto duplicate = std::find_if(std::next(first), methods.end(),
                                            [first](const std::unique_ptr<Method>& second) -> bool {
                                                return (*first)->id() == second->id();
                                            });

        if (duplicate != methods.end()) {
            throw std::invalid_argument{"method registry contains duplicate method id '" +
                                        std::string{(*first)->id()} + "'"};
        }
    }
}

} // namespace

MethodRegistry::MethodRegistry(std::vector<std::unique_ptr<Method>> methods)
    : methods_{std::move(methods)} {
    validate_methods(methods_);
    // List methods by descending priority, then ID.
    std::ranges::sort(
        methods_, [](const std::unique_ptr<Method>& first, const std::unique_ptr<Method>& second) {
            const auto first_priority = first->metadata().priority;
            const auto second_priority = second->metadata().priority;
            return first_priority != second_priority ? first_priority > second_priority
                                                     : first->id() < second->id();
        });
}

auto MethodRegistry::find(const std::string_view id) const noexcept -> const Method* {
    const auto iter =
        std::ranges::find_if(methods_, [id](const std::unique_ptr<Method>& method) -> bool {
            return method->id() == id;
        });

    if (iter == methods_.end()) {
        return nullptr;
    }

    return iter->get();
}

auto MethodRegistry::methods() const noexcept -> std::span<const std::unique_ptr<Method>> {
    return methods_;
}

auto method_registry() -> const MethodRegistry& {
    static const MethodRegistry registry{make_builtin_methods()};

    return registry;
}

} // namespace chargefw::methods
