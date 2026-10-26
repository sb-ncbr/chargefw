#include "support/test_methods.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/method_registry.h>

#include <memory>
#include <optional>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace methods = chargefw::methods;

namespace {} // namespace

TEST_CASE("method registry stores, finds, and rejects invalid methods",
          "[methods][method-registry]") {
    auto methods_for_registry = std::vector<std::unique_ptr<methods::Method>>{};
    methods_for_registry.push_back(std::make_unique<chargefw::test::StubMethod>("zeta"));
    methods_for_registry.push_back(std::make_unique<chargefw::test::StubMethod>("alpha"));
    methods_for_registry.push_back(
        std::make_unique<chargefw::test::StubMethod>("beta", methods::MethodRequirements{}, 5));

    const methods::MethodRegistry registry{std::move(methods_for_registry)};

    CHECK(registry.methods().size() == 3);
    CHECK(registry.find("alpha") != nullptr);
    CHECK(registry.find("zeta") != nullptr);
    CHECK(registry.find("missing") == nullptr);
    // Selection order: higher priority first, then ID.
    CHECK(registry.methods()[0]->id() == "beta");
    CHECK(registry.methods()[1]->id() == "alpha");
    CHECK(registry.methods()[2]->id() == "zeta");

    const auto null_method = [] {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::unique_ptr<methods::Method>{});
        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    };
    CHECK_THROWS_AS(null_method(), std::invalid_argument);

    const auto empty_id = [] {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::make_unique<chargefw::test::StubMethod>(""));
        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    };
    CHECK_THROWS_AS(empty_id(), std::invalid_argument);

    const auto duplicate_id = [] {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::make_unique<chargefw::test::StubMethod>("same"));
        invalid_methods.push_back(std::make_unique<chargefw::test::StubMethod>("same"));
        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    };
    CHECK_THROWS_AS(duplicate_id(), std::invalid_argument);
}
