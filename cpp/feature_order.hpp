#pragma once
#include <cstdint>
#include <utility>
#include <vector>

namespace stl2curves::detail {
// Mesh.nbrs in the reference is an integer set, subsequently walked as a list.
// Preserve CPython 3.12's insertion/probe/resize order for these nonnegative
// indices. Only insertion and iteration are needed; no Python runtime is used.
class FeatureNeighborOrder {
    std::vector<std::int64_t> slots = std::vector<std::int64_t>(8, -1);
    std::size_t used = 0;

    bool insert(unsigned key) {
        const std::uint64_t mask = slots.size() - 1;
        std::uint64_t at = key & mask, perturb = key;
        for (;;) {
            const unsigned count = at + 9 <= mask ? 10 : 1;
            for (unsigned offset = 0; offset < count; ++offset) {
                auto &slot = slots[at + offset];
                if (slot == key)
                    return false;
                if (slot < 0) {
                    slot = key;
                    return true;
                }
            }
            perturb >>= 5;
            at = (5 * at + 1 + perturb) & mask;
        }
    }

  public:
    void add(unsigned key) {
        if (!insert(key))
            return;
        ++used;
        if (used * 5 < (slots.size() - 1) * 3)
            return;
        const auto minimum = used * (used > 50000 ? 2 : 4);
        std::size_t size = 8;
        while (size <= minimum)
            size *= 2;
        auto previous = std::move(slots);
        slots.assign(size, -1);
        for (auto value : previous)
            if (value >= 0)
                insert(unsigned(value));
    }
    std::vector<unsigned> values() const {
        std::vector<unsigned> result;
        result.reserve(used);
        for (auto value : slots)
            if (value >= 0)
                result.push_back(unsigned(value));
        return result;
    }
};
} // namespace stl2curves::detail
