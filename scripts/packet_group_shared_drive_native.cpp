#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace {

constexpr int kStates = 65;

struct PointCapBuilder {
    int group_bits;
    int atoms;
    const double* local;
    double* maxima;
    std::array<std::array<double, kStates>, kStates> polynomials{};

    void recurse(int depth, int minimum_weight, int total_drive) {
        if (depth == atoms) {
            const int degree = depth * group_bits;
            for (int state_weight = 0; state_weight <= degree; ++state_weight) {
                double& destination =
                    maxima[state_weight * kStates + total_drive];
                destination = std::max(destination, polynomials[depth][state_weight]);
            }
            return;
        }

        const int degree = depth * group_bits;
        for (int weight = minimum_weight; weight <= group_bits; ++weight) {
            auto& next = polynomials[depth + 1];
            std::fill(next.begin(), next.end(), 0.0);
            const double* factor = local + weight * (group_bits + 1);
            for (int left = 0; left <= degree; ++left) {
                const double coefficient = polynomials[depth][left];
                if (coefficient == 0.0) {
                    continue;
                }
                for (int right = 0; right <= group_bits; ++right) {
                    next[left + right] += coefficient * factor[right];
                }
            }
            recurse(depth + 1, weight, total_drive + weight);
        }
    }
};

constexpr std::size_t source_offset(int incoming, int emitted, int source) {
    return (static_cast<std::size_t>(incoming) * kStates + emitted) * kStates +
           source;
}

constexpr std::size_t split_offset(int emitted, int state) {
    return static_cast<std::size_t>(emitted) * kStates + state;
}

}  // namespace

#if defined(_WIN32)
#define PACKET_GROUP_EXPORT extern "C" __declspec(dllexport)
#else
#define PACKET_GROUP_EXPORT extern "C"
#endif

// Enumerate nondecreasing drive-weight multisets with a shared-prefix
// convolution tree.  The Python reference rebuilds every polynomial from
// scratch for each leaf; sharing prefixes is especially valuable at g=8.
// This kernel is diagnostic binary64 only.  The outward certificate keeps its
// separately rounded Python implementation.
PACKET_GROUP_EXPORT int packet_group_point_caps(
    int group_bits,
    const double* local,
    double* maxima) {
    if (group_bits <= 0 || group_bits > 64 || 64 % group_bits != 0) {
        return 1;
    }
    std::fill(maxima, maxima + kStates * kStates, 0.0);
    PointCapBuilder builder{
        group_bits,
        64 / group_bits,
        local,
        maxima,
    };
    builder.polynomials[0][0] = 1.0;
    builder.recurse(0, 0, 0);
    return 0;
}

// The Python caller computes the exact NumPy argsort order.  Passing that
// order into this hot loop keeps tie handling identical to the reference
// implementation while moving the millions of scalar transport iterations
// out of the interpreter.
PACKET_GROUP_EXPORT int packet_group_shared_drive_apply(
    const double* source_caps,
    const double* source_units,
    const std::int32_t* source_counts,
    const double* split_caps,
    const double* pole_powers,
    const double* values,
    const std::int32_t* order,
    double* result) {
    for (int incoming = 0; incoming < kStates; ++incoming) {
        double total = 0.0;
        for (int emitted = 0; emitted < kStates; ++emitted) {
            const int count = source_counts[incoming * kStates + emitted];
            if (count == 0) {
                continue;
            }

            double source_mass = 0.0;
            for (int source = 0; source < count; ++source) {
                const std::size_t offset =
                    source_offset(incoming, emitted, source);
                source_mass += source_caps[offset] * source_units[offset];
            }

            int source_index = 0;
            int sink_index = 0;
            double source_remaining =
                source_units[source_offset(incoming, emitted, 0)];
            double sink_remaining =
                split_caps[split_offset(emitted, order[0])];
            double transported_value = 0.0;

            while (source_index < count) {
                while (sink_index < kStates && sink_remaining <= 0.0) {
                    ++sink_index;
                    if (sink_index < kStates) {
                        sink_remaining = split_caps[
                            split_offset(emitted, order[sink_index])];
                    }
                }
                if (sink_index == kStates) {
                    double residual_mass =
                        source_remaining *
                        source_caps[source_offset(
                            incoming, emitted, source_index)];
                    for (int source = source_index + 1; source < count;
                         ++source) {
                        const std::size_t offset =
                            source_offset(incoming, emitted, source);
                        residual_mass +=
                            source_caps[offset] * source_units[offset];
                    }
                    if (residual_mass >
                        1e-10 * std::max(1.0, source_mass)) {
                        return 1;
                    }
                    transported_value += residual_mass * values[order[0]];
                    source_index = count;
                    break;
                }

                const double take = std::min(source_remaining, sink_remaining);
                transported_value +=
                    take *
                    source_caps[source_offset(
                        incoming, emitted, source_index)] *
                    values[order[sink_index]];
                if (source_remaining <= sink_remaining) {
                    sink_remaining -= source_remaining;
                    ++source_index;
                    if (source_index < count) {
                        source_remaining = source_units[source_offset(
                            incoming, emitted, source_index)];
                    }
                } else {
                    source_remaining -= sink_remaining;
                    ++sink_index;
                    if (sink_index < kStates) {
                        sink_remaining = split_caps[
                            split_offset(emitted, order[sink_index])];
                    }
                }
            }
            total += pole_powers[emitted] * transported_value;
        }
        result[incoming] = total;
    }
    return 0;
}
