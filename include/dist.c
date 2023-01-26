#include "dist.h"

#include "util.h"

#include <flint/fmpz_poly.h>
#include <stdbool.h>

//#region global definitions for Ziggurat algorithm
#define SHR3 (var_j = jsr, jsr ^= (jsr << 13), jsr ^= (jsr >> 17), jsr ^= (jsr << 5),var_j + jsr)
#define UNI (0.5 + (signed) SHR3 * 0.2328306e-9)
#define RNOR (var_h = SHR3, var_i = var_h & 127, (abs(var_h)<k[var_i])? var_h*w[var_i] : ziggurat_fallback())

static unsigned long var_i, var_j, jsr;
static unsigned long k[128];
static int var_h;
static double w[128], f[128];
//#endregion Required for Ziggurat algorithm

/// Get a uniformly random number in the range (0,1)
static inline __attribute__((always_inline)) double get_uniformly_random() {
    //Get 32 bits of randomness
    //32 bit is enough as the used std_deviation is 8
    //This is not a security issue because the security does only depend on gaussian distribution, not on bit size
    //i.e. this does not become a problem if the Gaussian value does not exceed 2^32

    static double max = (double) ((1L << 32) - 1);
    double r;

    mpz_t rand;
    mpz_init2(rand, 32);

    get_random(rand, 32);

    r = (mpz_get_d(rand) / max);
    mpz_clear(rand);

    return r;
}

double dist_gauss_box_muller(const double std_deviation) {
    // z0 = std_dev * sqrt(-2 log(u1)) cos (2 pi u2)
    // z1 = std_dev * sqrt(-2 log(u1)) sin (2 pi u2)
    // r = sqrt(-2 log(u1))
    // phi = 2 pi u2
    // -> z0 = std_dev * r * cos(phi)
    // -> z1 = std_dev * r * sin(phi)

    static bool cached = false;
    static double z1;

    if (cached == false) {
        //No value cached ,generate random values
        double r, phi;

        r = sqrt(-2.0 * log(get_uniformly_random()));
        phi = 2.0 * M_PI * (get_uniformly_random());

        z1 = std_deviation * r * sin(phi);
        cached = true;

        return (std_deviation * r * cos(phi));
    }

    cached = false;
    return z1;
}

double dist_gauss_polar(const double std_deviation) {
    //r^2=y_1^2 + y_2^2 < 1
    //x_1 = sqrt(-2 * log(r^2)/r^2) * y_1
    //x_2 = sqrt(-2 * log(r^2)/r^2) * y_2

    static bool cached = false;
    static double x2;

    if (cached == false) {
        double y1, y2, r_2, t;

        do {
            y1 = get_uniformly_random() * 2.0 - 1.0;
            y2 = get_uniformly_random() * 2.0 - 1.0;
            r_2 = y1 * y1 + y2 * y2;
        } while (r_2 >= 1.0 || r_2 == 0);

        t = sqrt(-2.0 * log(r_2) / r_2) * std_deviation;

        x2 = y2 * t;

        cached = true;

        return y1 * t;;
    }

    cached = false;
    return x2;
}

/// Fallback algorithm for RNOR #define
static double ziggurat_fallback() {
    const double r = 3.442620f;
    static double x, y;
    for (;;) {
        x = (double) var_h * w[var_i];

        if (var_i == 0) {
            do {
                x = -log(UNI) * 0.2904764;
                y = -log((UNI));
            } while (y + y < x * x);

            if (var_h > 0) {
                return r + x;
            } else {
                return -r - x;
            }
        }

        if (f[var_i] + UNI * (f[var_i - 1] - f[var_i]) < exp(-.5 * x * x)) {
            return x;
        }

        var_h = SHR3;
        var_i = var_h & 127;

        if (abs(var_h) < k[var_i])
            return ((float) var_h * w[var_i]);
    }
}

/// Initialize the ziggurat sampler, precompute tables
static void zigset() {
    const double m = 2147483648.0;

    double d = 3.442619855899;
    double t = d;
    double v = 9.91256303526217e-3;
    double q;

    int i;

    mpz_t rand;
    mpz_init(rand);
    get_random(rand, 32);
    jsr = mpz_get_ui(rand);
    mpz_clear(rand);

    /* Tables for RNOR: */ q = v / exp(-.5 * d * d);
    k[0] = (unsigned long) ((d / q) * m);
    k[1] = 0;
    w[0] = (float) (q / m);
    w[127] = (float) (d / m);
    f[0] = 1.0f;
    f[127] = exp(-.5 * d * d);
    for (i = 126; i >= 1; i--) {
        d = sqrt(-2. * log(v / d + exp(-.5 * d * d)));
        k[i + 1] = (unsigned long) ((d / t) * m);
        t = d;
        f[i] = exp(-0.5 * d * d);
        w[i] = d / m;
    }
}

double dist_gauss_ziggurat(const double std_deviation) {
    static bool initialized = false;

    if (initialized == false){
        zigset();
        initialized = true;
    }

    return RNOR * std_deviation;
}