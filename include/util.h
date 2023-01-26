#ifndef CUSTOM_UTIL_H
#define CUSTOM_UTIL_H

#define INT_SIZE 4          // sizeof(int)
#define INT_BIT_SIZE 32     // sizeof(int) * 8

#include <flint/fmpz_poly.h>

struct settings {
    signed long n;     // Degree of polynomials
    fmpz_t q;               // Large prime
    unsigned long qBits;        // Bits of q; effectively lg2(q)
    unsigned long t;            // Message space
    signed int b;           // Encoding base
    unsigned long D;            // Ciphertext max_len/Maximum degree of homomorphism
    double std_dev;         // Standard deviation of gaussian distribution
    double greater_std_dev; // Greater standard deviation of gaussian distribution
};

/// Fetch count * 32 random bits
/// IMPORTANT This function does not check array max_len
/// and might write to arbitrary memory if count > sizeof(data)
/// @param[out] data Empty array of size >= count
/// @param[in] count Amount of 4 byte blocks to fetch
void urandom(unsigned int data[], unsigned long count);

/// Get n bit random values
/// @param[out] random Empty Arbitrary precision integer to take random data
/// @param[in] bits Amount of bits
void get_random(mpz_t random, unsigned long bits);

/// Generate a random large prime of a certain bit-size
/// @param[out] prime Empty FLINT Arbitrary precision integer to take random prime
/// @param[in] bits Bit-size
void generate_prime(fmpz_t prime, unsigned long bits);

/// Generate a random large prime of a certain bit-size where 1 = q mod 2n
/// @param[out] prime Empty FLINT Arbitrary precision integer to take random prime
/// @param[in] bits Minimum bit-size
/// @param[in] n Polynomial degree n
void generate_prime_congruent_mod_2n(fmpz_t prime, unsigned long bits, unsigned long n);

/// Fill settings with parameters
/// @param[out] settings Empty settings
/// @param[in] n_power Power of n
/// @param[in] q Modulus q
/// @param[in] t Message space
/// @param[in] b Encoding base
/// @param[in] D Maximum ciphertext length / maximum homomorphic depth minus 2
void settings_init(struct settings *settings, unsigned long n_power, fmpz_t q, unsigned long t, signed int b, unsigned long D);

/// Check settings for errors
/// @param[in] settings Settings
/// @return 0 if struct is ok, != 0 otherwise
int settings_check(struct settings settings);

/// Print settings
/// @param[in] settings Settings
void settings_print(struct settings settings);

/// Generate standard deviation
/// @param[in] n Polynomial degree n
/// @return Standard deviation
double gen_std_deviation(signed long n);

/// Generate greater standard deviation
/// @param[in] std_dev Standard deviation
/// @param[in] n Polynomial degree n
/// @return Greater standard deviation
double gen_greater_std_deviation(double std_dev, signed long n);

#endif //CUSTOM_UTIL_H
