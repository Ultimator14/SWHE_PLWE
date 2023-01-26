#include "wrapper.h"

#include "asym.h"
#include "encoding.h"
#include "plwe_poly.h"
#include "util.h"
#include "message.h"

#include <flint/fmpz_poly.h>

void settings_init_gen_prime(struct settings *settings, unsigned long n_power, unsigned long qBits, unsigned long t, signed int b, unsigned long D) {
    fmpz_t q;
    fmpz_init(q);

    generate_prime(q, qBits);
    settings_init(settings, n_power, q, t, b, D);

    fmpz_clear(q);
}

void settings_init_gen_prime_congruent_mod_2n(struct settings *settings, unsigned long n_power, unsigned long qBits, unsigned long t, signed int b, unsigned long D) {
    fmpz_t q;
    fmpz_init(q);

    generate_prime_congruent_mod_2n(q, qBits, settings->n);
    settings_init(settings, n_power, q, t, b, D);

    fmpz_clear(q);
}

void encode_encrypt(struct message *output, signed long input, const struct settings *settings, const struct key *key){
    mpz_t in;
    mpz_init_set_si(in, input);

    struct plwe_poly poly;
    plwe_poly_init(&poly, settings->q, settings->n);

    encode(&poly, in, settings->b);
    encrypt(output, &poly, key);

    plwe_poly_clear(&poly);
    mpz_clear(in);
}

signed int decrypt_decode(struct message *input, const struct settings *settings, const struct key *key){
    mpz_t out;
    mpz_init(out);

    struct plwe_poly poly;
    plwe_poly_init(&poly, settings->q, settings->n);

    decrypt(&poly, input, key);
    decode(out, &poly, settings->b);

    signed int ret = mpz_get_si(out);
    mpz_clear(out);

    return ret;

}

void encode_eval_add_plain(struct message *output, struct message *message, const signed long plain, const struct settings *settings) {
    mpz_t in;
    mpz_init_set_si(in, plain);

    struct plwe_poly poly;
    plwe_poly_init(&poly, settings->q, settings->n);
    encode(&poly, in, settings->b);

    eval_add_plain(output, message, &poly);

    mpz_clear(in);
}

void encode_eval_mul_plain(struct message *output, struct message *message, const signed long plain, const struct settings *settings) {
    mpz_t in;
    mpz_init_set_si(in, plain);

    struct plwe_poly poly;
    plwe_poly_init(&poly, settings->q, settings->n);
    encode(&poly, in, settings->b);

    eval_mul_plain(output, message, &poly);

    mpz_clear(in);
}
