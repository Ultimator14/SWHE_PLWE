#ifndef CUSTOM_ENCODING_H
#define CUSTOM_ENCODING_H

#include <flint/fmpz_poly.h>

//Forward declarations
struct plwe_poly;   /// defined in plwe_poly.h

/// Encode an arbitrary sized integer to a polynomial
/// @param[out] output Polynomial
/// @param[in] input Integer
/// @param[in] b Base
void encode(struct plwe_poly *output, const mpz_t input, signed int b);

/// Decode a plwe_poly to a mpz_t

/// Decode a polynomial to an arbitrary sized integer
/// @param[out] output Integer
/// @param[in] input Polynomial
/// @param[in] b Base
void decode(mpz_t output, const struct plwe_poly *input, signed int b);

#endif //CUSTOM_ENCODING_H
