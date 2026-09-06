// Log-domain pair-support DP. Binary64 diagnostic, never an outward bound.
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace {
constexpr double neg_inf = -std::numeric_limits<double>::infinity();
inline double add(double a, double b) {
    if (a == neg_inf) return b;
    if (b == neg_inf) return a;
    const double hi = a > b ? a : b;
    const double lo = a > b ? b : a;
    return hi + std::log1p(std::exp(lo - hi));
}
inline void product(const double* v, const double* m, double* out) {
    out[0] = add(add(v[0]+m[0], v[1]+m[3]), v[2]+m[6]);
    out[1] = add(add(v[0]+m[1], v[1]+m[4]), v[2]+m[7]);
    out[2] = add(add(v[0]+m[2], v[1]+m[5]), v[2]+m[8]);
}
inline void accumulate(const double* v, const double* m, double* out) {
    double term[3];
    product(v, m, term);
    out[0] = add(out[0], term[0]);
    out[1] = add(out[1], term[1]);
    out[2] = add(out[2], term[2]);
}
}

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

extern "C" EXPORT int pair_log_coefficients(const double* region, int length, double* output) noexcept {
    if (!region || !output || length < 1 || length > 2048) return 1;
    try {
        const std::size_t stride = static_cast<std::size_t>(length) + 1;
        std::vector<double> first(stride*stride*3, neg_inf);
        std::vector<double> second(stride*stride*3, neg_inf);
        double* current = first.data();
        double* next = second.data();
        current[0] = 0.;
        for (int n=0; n<length; ++n) {
            // Symmetry exchanges the two independently shuffled outer rows.
            // Only the lower triangle is evaluated; mirror each three-state cell.
            for (int a=0; a<=n+1; ++a) {
                for (int b=0; b<=a; ++b) {
                    double value[3] = {neg_inf, neg_inf, neg_inf};
                    if (a<=n && b<=n)
                        product(current+(a*stride+b)*3, region, value);
                    if (a && b<=n)
                        accumulate(current+((a-1)*stride+b)*3, region+9, value);
                    if (b && a<=n)
                        accumulate(current+(a*stride+b-1)*3, region+9, value);
                    if (a && b)
                        accumulate(current+((a-1)*stride+b-1)*3, region+18, value);
                    double* left = next+(a*stride+b)*3;
                    double* right = next+(b*stride+a)*3;
                    left[0]=right[0]=value[0];
                    left[1]=right[1]=value[1];
                    left[2]=right[2]=value[2];
                }
            }
            double* swap=current; current=next; next=swap;
        }
        for (std::size_t i=0; i<stride*stride; ++i)
            output[i]=add(add(current[3*i], current[3*i+1]), current[3*i+2]);
        return 0;
    } catch (...) { return 2; }
}
