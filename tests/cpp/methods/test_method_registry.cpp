#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>

#include <cstddef>
#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("method registry finds built-in methods in selection order",
          "[methods][method-registry]") {
    const auto& registry = methods::method_registry();

    REQUIRE_FALSE(registry.methods().empty());
    CHECK(registry.find("eem") != nullptr);
    CHECK(registry.find("missing") == nullptr);

    for (const auto& method : registry.methods()) {
        CAPTURE(method->id());
        CHECK(registry.find(method->id()) == method.get());
    }

    // Selection order: higher priority first, then ID.
    for (std::size_t index = 1; index < registry.methods().size(); ++index) {
        const auto& previous = *registry.methods()[index - 1];
        const auto& current = *registry.methods()[index];
        CAPTURE(previous.id(), current.id());
        CHECK((previous.metadata().priority > current.metadata().priority ||
               (previous.metadata().priority == current.metadata().priority &&
                previous.id() < current.id())));
    }
}
