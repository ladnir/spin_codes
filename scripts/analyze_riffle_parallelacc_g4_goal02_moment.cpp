#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t kPacketCount = 524352;
constexpr std::uint32_t kDistance = 188766;
constexpr long double kLn2 = 0.693147180559945309417232121458176568L;

struct Eval {
    long double objective;
    std::array<long double, 2> gradient;
};

struct ProfileDp {
    std::string id;
    std::array<std::uint8_t, 15> counts{};
    std::vector<std::uint8_t> values;
    std::vector<std::uint8_t> maxima;
    std::vector<std::uint32_t> strides;
    std::vector<std::vector<std::uint8_t>> remaining;
    std::vector<std::uint8_t> xor_state;
    std::vector<long double> mass;
    std::vector<long double> derivative_a;
    std::vector<long double> derivative_b;
    std::uint32_t state_count = 0;
    std::uint32_t support = 0;
    std::uint64_t evaluations = 0;

    explicit ProfileDp(std::string profile_id, const std::array<std::uint8_t, 15>& profile_counts)
        : id(std::move(profile_id)), counts(profile_counts) {
        std::uint64_t product = 1;
        for (std::uint32_t value = 1; value <= 15; ++value) {
            const auto count = counts[value - 1];
            support += count;
            if (count == 0) {
                continue;
            }
            values.push_back(static_cast<std::uint8_t>(value));
            maxima.push_back(count);
            strides.push_back(static_cast<std::uint32_t>(product));
            product *= static_cast<std::uint64_t>(count) + 1;
        }
        if (support == 0 || product > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("invalid profile state space");
        }
        state_count = static_cast<std::uint32_t>(product);
        remaining.resize(values.size(), std::vector<std::uint8_t>(state_count));
        xor_state.resize(state_count);
        mass.resize(state_count);
        derivative_a.resize(state_count);
        derivative_b.resize(state_count);

        for (std::uint32_t index = 0; index < state_count; ++index) {
            std::uint8_t state = 0;
            for (std::size_t slot = 0; slot < values.size(); ++slot) {
                const auto radix = static_cast<std::uint32_t>(maxima[slot]) + 1;
                const auto used = static_cast<std::uint8_t>((index / strides[slot]) % radix);
                remaining[slot][index] = static_cast<std::uint8_t>(maxima[slot] - used);
                if (used & 1U) {
                    state ^= values[slot];
                }
            }
            xor_state[index] = state;
        }
    }

    Eval evaluate(long double u, long double v) {
        ++evaluations;
        const long double scale = static_cast<long double>(support) / kPacketCount;
        const long double a = scale * std::exp(u);
        const long double b = scale * std::exp(v);
        const long double tau = std::exp(-b);

        std::array<long double, 5> factor{};
        std::array<long double, 5> dlog_da{};
        std::array<long double, 5> dlog_db{};
        for (std::uint32_t weight = 0; weight <= 4; ++weight) {
            const long double x = std::exp(-b - a * weight);
            const long double denominator = -std::expm1(-b - a * weight);
            factor[weight] = std::exp(-a * weight) / denominator;
            dlog_da[weight] = -static_cast<long double>(weight) / denominator;
            dlog_db[weight] = -x / denominator;
        }

        std::fill(mass.begin(), mass.end(), 0.0L);
        std::fill(derivative_a.begin(), derivative_a.end(), 0.0L);
        std::fill(derivative_b.begin(), derivative_b.end(), 0.0L);
        mass[0] = 1.0L;

        for (std::uint32_t index = 0; index < state_count; ++index) {
            const long double current_mass = mass[index];
            if (current_mass == 0.0L) {
                continue;
            }
            const auto current_state = xor_state[index];
            for (std::size_t slot = 0; slot < values.size(); ++slot) {
                const auto multiplicity = remaining[slot][index];
                if (multiplicity == 0) {
                    continue;
                }
                const auto next = index + strides[slot];
                const auto next_state = static_cast<std::uint8_t>(current_state ^ values[slot]);
                const auto weight = static_cast<std::uint32_t>(std::popcount(next_state));
                const long double transition = static_cast<long double>(multiplicity) * factor[weight];
                mass[next] += current_mass * transition;
                derivative_a[next] +=
                    (derivative_a[index] + current_mass * dlog_da[weight]) * transition;
                derivative_b[next] +=
                    (derivative_b[index] + current_mass * dlog_db[weight]) * transition;
            }
        }

        const auto final_index = state_count - 1;
        const long double final_mass = mass[final_index];
        if (!(final_mass > 0.0L) || !std::isfinite(final_mass)) {
            return {std::numeric_limits<long double>::infinity(), {0.0L, 0.0L}};
        }

        const long double log_choose =
            std::lgammal(static_cast<long double>(kPacketCount) + 1.0L) -
            std::lgammal(static_cast<long double>(support) + 1.0L) -
            std::lgammal(static_cast<long double>(kPacketCount - support) + 1.0L);
        const long double objective =
            static_cast<long double>(kDistance) * a +
            static_cast<long double>(kPacketCount - support) * b - log_choose -
            std::log1pl(-tau) - std::lgammal(static_cast<long double>(support) + 1.0L) +
            std::log(final_mass);

        const long double gradient_a =
            static_cast<long double>(kDistance) + derivative_a[final_index] / final_mass;
        const long double gradient_b =
            static_cast<long double>(kPacketCount - support) - tau / (1.0L - tau) +
            derivative_b[final_index] / final_mass;
        return {objective, {a * gradient_a, b * gradient_b}};
    }
};

struct Optimum {
    long double objective;
    long double u;
    long double v;
    std::array<long double, 2> gradient;
};

Optimum optimize(ProfileDp& profile) {
    Optimum best{std::numeric_limits<long double>::infinity(), 0.0L, 0.0L, {0.0L, 0.0L}};
    for (const long double u : {-2.0L, 0.0L, 2.0L}) {
        for (const long double v : {-2.0L, 0.0L, 2.0L}) {
            const auto evaluated = profile.evaluate(u, v);
            if (evaluated.objective < best.objective) {
                best = {evaluated.objective, u, v, evaluated.gradient};
            }
        }
    }

    long double h00 = 1.0L;
    long double h01 = 0.0L;
    long double h11 = 1.0L;
    for (std::uint32_t iteration = 0; iteration < 40; ++iteration) {
        const long double gradient_norm = std::hypotl(best.gradient[0], best.gradient[1]);
        if (gradient_norm < 1e-10L) {
            break;
        }

        long double p0 = -(h00 * best.gradient[0] + h01 * best.gradient[1]);
        long double p1 = -(h01 * best.gradient[0] + h11 * best.gradient[1]);
        long double directional = best.gradient[0] * p0 + best.gradient[1] * p1;
        if (!(directional < 0.0L)) {
            p0 = -best.gradient[0];
            p1 = -best.gradient[1];
            directional = -(best.gradient[0] * best.gradient[0] + best.gradient[1] * best.gradient[1]);
            h00 = h11 = 1.0L;
            h01 = 0.0L;
        }
        const long double direction_norm = std::hypotl(p0, p1);
        if (direction_norm > 2.0L) {
            p0 *= 2.0L / direction_norm;
            p1 *= 2.0L / direction_norm;
            directional = best.gradient[0] * p0 + best.gradient[1] * p1;
        }

        long double step = 1.0L;
        Optimum candidate = best;
        bool accepted = false;
        for (std::uint32_t line = 0; line < 16; ++line) {
            const long double next_u = best.u + step * p0;
            const long double next_v = best.v + step * p1;
            const auto evaluated = profile.evaluate(next_u, next_v);
            if (evaluated.objective <= best.objective + 1e-4L * step * directional) {
                candidate = {evaluated.objective, next_u, next_v, evaluated.gradient};
                accepted = true;
                break;
            }
            step *= 0.5L;
        }
        if (!accepted) {
            break;
        }

        const long double s0 = candidate.u - best.u;
        const long double s1 = candidate.v - best.v;
        const long double y0 = candidate.gradient[0] - best.gradient[0];
        const long double y1 = candidate.gradient[1] - best.gradient[1];
        const long double ys = y0 * s0 + y1 * s1;
        if (ys > 1e-14L) {
            const long double rho = 1.0L / ys;
            const long double hy0 = h00 * y0 + h01 * y1;
            const long double hy1 = h01 * y0 + h11 * y1;
            const long double yhy = y0 * hy0 + y1 * hy1;
            const long double coefficient = (1.0L + rho * yhy) * rho;
            h00 += coefficient * s0 * s0 - rho * (s0 * hy0 + hy0 * s0);
            h01 += coefficient * s0 * s1 - rho * (s0 * hy1 + hy0 * s1);
            h11 += coefficient * s1 * s1 - rho * (s1 * hy1 + hy1 * s1);
        } else {
            h00 = h11 = 1.0L;
            h01 = 0.0L;
        }
        best = candidate;
    }
    return best;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::string id;
    while (std::cin >> id) {
        std::array<std::uint8_t, 15> counts{};
        for (auto& count : counts) {
            std::uint32_t parsed = 0;
            if (!(std::cin >> parsed) || parsed > 255) {
                std::cerr << "invalid profile line\n";
                return 2;
            }
            count = static_cast<std::uint8_t>(parsed);
        }
        ProfileDp profile(id, counts);
        const auto optimum = optimize(profile);
        const long double scale = static_cast<long double>(profile.support) / kPacketCount;
        const long double a = scale * std::exp(optimum.u);
        const long double b = scale * std::exp(optimum.v);
        std::cout << id << '\t' << profile.support << '\t' << profile.state_count << '\t'
                  << profile.evaluations << '\t' << std::setprecision(18)
                  << optimum.objective / kLn2 << '\t' << a << '\t' << b << '\t'
                  << optimum.gradient[0] << '\t' << optimum.gradient[1] << '\n';
    }
    return 0;
}
