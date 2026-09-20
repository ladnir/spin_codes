// Four-state uniform-refresh Q1. Nearest binary64, not outward arithmetic.
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
    out[0] = add(add(v[0]+m[0], v[1]+m[4]), add(v[2]+m[8], v[3]+m[12]));
    out[1] = add(add(v[0]+m[1], v[1]+m[5]), add(v[2]+m[9], v[3]+m[13]));
    out[2] = add(add(v[0]+m[2], v[1]+m[6]), add(v[2]+m[10], v[3]+m[14]));
    out[3] = add(add(v[0]+m[3], v[1]+m[7]), add(v[2]+m[11], v[3]+m[15]));
}
}

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

extern "C" EXPORT int refresh_log_coefficients(const double* zero, const double* one,
    int batch, int length, double* output) noexcept {
    if (!zero || !one || !output || batch<1 || length<1 || length>4096) return 1;
    try {
        const std::size_t stride = static_cast<std::size_t>(length)+1;
        std::vector<double> first(stride*4);
        std::vector<double> second(stride*4);
        for (int z=0; z<batch; ++z) {
            double* current=first.data();
            double* next=second.data();
            current[0]=0.; current[1]=current[2]=current[3]=neg_inf;
            const double* rz=zero+z*16;
            const double* ra=one+z*16;
            for (int n=0; n<length; ++n) {
                product(current,rz,next);
                for (int w=1; w<=n; ++w) {
                    double term[4];
                    product(current+w*4,rz,next+w*4);
                    product(current+(w-1)*4,ra,term);
                    next[w*4]   = add(next[w*4],term[0]);
                    next[w*4+1] = add(next[w*4+1],term[1]);
                    next[w*4+2] = add(next[w*4+2],term[2]);
                    next[w*4+3] = add(next[w*4+3],term[3]);
                }
                product(current+n*4,ra,next+(n+1)*4);
                double* swap=current; current=next; next=swap;
            }
            for (std::size_t w=0; w<stride; ++w)
                output[z*stride+w]=add(add(current[w*4],current[w*4+1]),add(current[w*4+2],current[w*4+3]));
        }
        return 0;
    } catch (...) { return 2; }
}
