// Batched Q1 coefficient recurrence. Nearest binary64, not outward arithmetic.
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
    return hi + std::log1p(std::exp(lo-hi));
}
inline void product(const double* v, const double* m, double* out) {
    out[0] = add(add(v[0]+m[0], v[1]+m[3]), v[2]+m[6]);
    out[1] = add(add(v[0]+m[1], v[1]+m[4]), v[2]+m[7]);
    out[2] = add(add(v[0]+m[2], v[1]+m[5]), v[2]+m[8]);
}
}

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

extern "C" EXPORT int q1_log_coefficients(const double* zero, const double* one,
    int batch, int length, double* output) noexcept {
    if (!zero || !one || !output || batch<1 || length<1 || length>4096) return 1;
    try {
        const std::size_t stride = static_cast<std::size_t>(length)+1;
        std::vector<double> first(stride*3);
        std::vector<double> second(stride*3);
        for (int z=0; z<batch; ++z) {
            double* current=first.data();
            double* next=second.data();
            current[0]=0.; current[1]=current[2]=neg_inf;
            const double* rz=zero+z*9;
            const double* ra=one+z*9;
            for (int n=0; n<length; ++n) {
                product(current,rz,next);
                for (int w=1; w<=n; ++w) {
                    double term[3];
                    product(current+w*3,rz,next+w*3);
                    product(current+(w-1)*3,ra,term);
                    next[w*3]   = add(next[w*3],term[0]);
                    next[w*3+1] = add(next[w*3+1],term[1]);
                    next[w*3+2] = add(next[w*3+2],term[2]);
                }
                product(current+n*3,ra,next+(n+1)*3);
                double* swap=current; current=next; next=swap;
            }
            for (std::size_t w=0; w<stride; ++w)
                output[z*stride+w]=add(add(current[w*3],current[w*3+1]),current[w*3+2]);
        }
        return 0;
    } catch (...) { return 2; }
}
