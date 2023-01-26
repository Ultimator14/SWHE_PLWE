#ifndef CUSTOM_PLWE_POLY_H
#define CUSTOM_PLWE_POLY_H

#include <flint/fmpz_poly.h>

struct plwe_poly {
    signed long n;
    fmpz_t mod;
    fmpz_poly_t fmod;
    fmpz_poly_t poly;
};

/// Initialize polynomial
/// @param[out] poly Empty polynomial
/// @param[in] q Coefficient modulus q
/// @param[in] n Polynomial degree n used for f(x)=x^n + 1
void plwe_poly_init(struct plwe_poly *poly, const fmpz_t q, signed long n);

/// Clear polynomial
/// @param[in,out] poly Polynomial
void plwe_poly_clear(struct plwe_poly *poly);

/// Reduce polynomial coefficients mod t
/// @param[in,out] poly Polynomial
/// @param[in] t Plaintext modulus t
void plwe_poly_mod_t(struct plwe_poly *poly, unsigned long t);

/// Reduce polynomial mod f(x) and coefficients mod q
/// @param[in,out] poly Polynomial
void plwe_poly_pmod(struct plwe_poly *poly);

/// Copy a polynomial (deep copy)
/// @param[out] out Target polynomial
/// @param[in] in Source polynomial
extern void plwe_poly_set(struct plwe_poly *out, const struct plwe_poly *in);

/// Add two polynomials
/// @param[out] result Result of the addition
/// @param[in] poly1 Polynomial 1
/// @param[in] poly2 Polynomial 2
extern void plwe_poly_add(struct plwe_poly *result, const struct plwe_poly *poly1, const struct plwe_poly *poly2);

/// Multiply two polynomials
/// @param[out] result Result of the multiplication
/// @param[in] poly1 Polynomial 1
/// @param[in] poly2 Polynomial 2
extern void plwe_poly_mul(struct plwe_poly *result, const struct plwe_poly *poly1, const struct plwe_poly *poly2);

/// Multiply a polynomial with an unsigned long integer
/// @param[out] result Result of the computation
/// @param[in] poly Polynomial
/// @param[in] scalar Scalar
extern void plwe_poly_scalar_mul_ui(struct plwe_poly *result, const struct plwe_poly *poly, unsigned long scalar);

/// Multiply a polynomial with a signed long integer
/// @param[out] result Result of the computation
/// @param[in] poly Polynomial
/// @param[in] scalar Scalar
extern void plwe_poly_scalar_mul_si(struct plwe_poly *result, const struct plwe_poly *poly, signed long scalar);

/// Multiply a polynomial with an arbitrary sized integer
/// @param[out] result Result of the computation
/// @param[in] poly Polynomial
/// @param[in] scalar Scalar
extern void plwe_poly_scalar_mul_mpz(struct plwe_poly *result, const struct plwe_poly *poly, mpz_t scalar);

/// Print a polynomial
/// @param[in] poly Polynomial
void plwe_poly_print(const struct plwe_poly *poly);

/// Generate a polynomial with uniformly distributed (random) coefficients
/// @param[out] poly Polynomial
/// @param[in] qBits Maximum bit-size of coefficients
void rand_poly_uniform(struct plwe_poly *poly, unsigned long qBits);

/// Generate a polynomial with gaussian distributed coefficients
/// @param[out] poly Polynomial
/// @param[in] std_dev Standard deviation of the gaussian distribution
void rand_poly_gauss(struct plwe_poly *poly, double std_dev);

#endif //CUSTOM_PLWE_POLY_H
