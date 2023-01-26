#include "util.h"

#ifdef LIB_SODIUM
#include <sodium.h>         // For random numbers
#endif
#ifndef LIB_SODIUM
//#include <stdio.h>          // For reading urandom
#endif

/// Compute the amount of blocks required for certain bit- and blocksize
/// @param[in] bits Amount of bits
/// @param[in] blocksize Blocksize (usually 4 bytes)
/// @return Amount of blocks
static inline __attribute__((always_inline)) unsigned long bits_to_blocks(unsigned long bits, unsigned long blocksize){
    return (bits + (blocksize * 8 - 1)) >> 5;  // bits to blocks of blocksize bytes (rounded up)
}

#ifdef LIB_SODIUM
void urandom(unsigned int data[], unsigned long count){
    randombytes_buf(data, count * INT_SIZE);
}
#endif
#ifndef LIB_SODIUM
void urandom(unsigned int data[], unsigned long count)
{
    FILE *fp;
    fp = fopen("/dev/urandom", "r");
    fread(data, INT_SIZE, count, fp);  //Read count * 4 Byte
    fclose(fp);
}
#endif

void get_random(mpz_t random, unsigned long bits){
    unsigned long block_count = bits_to_blocks(bits, 4);  // use integer array (4 bytes/block)
    unsigned int data[block_count];

    // Always fill whole array for simplicity, not only bit size
    urandom(data, block_count);

    mpz_import(random, block_count, 1, INT_SIZE, 1, 0, data);  //4*4 byte, most significant word first

    //Reuse the data variable and create bitmask for prime
    data[0] = (1UL << (((bits - 1) % INT_BIT_SIZE) + 1)) -1;
    for (int i = 1; i < block_count; i++) {
        data[i] = UINT_MAX;
    }

    mpz_t bitmask;
    mpz_init2(bitmask, bits);
    mpz_import(bitmask, block_count, 1, sizeof(unsigned int), 0, 0, data);

    mpz_and(random, random, bitmask);

    mpz_clear(bitmask);
}

void generate_prime(fmpz_t prime, unsigned long bits){
    mpz_t data;
    mpz_init2(data, bits);

    do {
        get_random(data, bits);
        mpz_nextprime(data, data);
    } while (mpz_sizeinbase(data, 2) != bits);  // repeat until used prime bits match bitsize

    fmpz_set_mpz(prime, data);
}

void generate_prime_congruent_mod_2n(fmpz_t prime, unsigned long bits, unsigned long n){
    generate_prime(prime, bits);

    mpz_t one, two_n, p;
    mpz_init_set_ui(one, 1);
    mpz_init_set_ui(two_n, n << 1);
    mpz_init(p);

    fmpz_get_mpz(p, prime);

    while (mpz_congruent_p(one, p, two_n) == 0) {
        mpz_nextprime(p, p);
    }

    fmpz_set_mpz(prime, p);

    mpz_clear(one);
    mpz_clear(two_n);
    mpz_clear(p);
}

void settings_init(struct settings *settings, unsigned long n_power, fmpz_t q, unsigned long t, signed int b, unsigned long D) {
    settings->n = 1 << n_power;  // n must be a power of 2
    fmpz_set(settings->q, q);
    settings->qBits = fmpz_sizeinbase(q, 2);
    settings->t = t;  // t must be co-prime to q (always true for all t<q since q is prime)
    settings->b = b;
    settings->D = D;

    settings->std_dev = gen_std_deviation(settings->n);
    settings->greater_std_dev = gen_greater_std_deviation((double) settings->std_dev, settings->n);
}

int settings_check(const struct settings settings)
{
    //Check n
    if ((settings.n == 0) || ((settings.n != 1) && settings.n & (settings.n - 1))) {
        return 1;
    }

    //Check qBits
    if (fmpz_sizeinbase(settings.q, 2) != settings.qBits) {
        return 2;
    }

    //Check q/qBits
    mpz_t pprime;
    mpz_init2(pprime,settings.qBits);
    fmpz_get_mpz(pprime, settings.q);

    if (! mpz_probab_prime_p(pprime, 500)) {
        return 3;
    }

    mpz_clear(pprime);

    //Check t
    if (fmpz_cmp_ui(settings.q, settings.t) < 0) {
        return 4;
    }

    //Check b
    if (settings.b < 2) {
        return 5;
    }

    return 0;
}

void settings_print(const struct settings settings) {
    printf("Settings\n"
           "----------------------------------\n");

    printf("n: %ld\n", settings.n);
    printf("q: ");
    fmpz_print(settings.q);
    printf("\n");
    printf("lg2(q): %ld\n", settings.qBits);
    printf("t: %ld\n", settings.t);
    printf("b: %d\n", settings.b);
    printf("D: %ld\n", settings.D);

    printf("----------------------------------\n");
}

double gen_std_deviation(signed long n){
    //Usually 2 ^(log n) i.e. 8 for n=1024, but 8 seems to be ok for security overall
    return 8;
}

double gen_greater_std_deviation(double std_dev, signed long n) {
    //Usually 2 ^(log n) * std_deviation, but we use lg2 -> n * std_deviation
    return std_dev * (double) n;
}
