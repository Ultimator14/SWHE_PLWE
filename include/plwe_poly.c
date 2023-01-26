#include "plwe_poly.h"

#include "dist.h"
#include "util.h"

void plwe_poly_init(struct plwe_poly *poly, const fmpz_t q, const signed long n) {
    // q = coefficient modulo
    // n = polynomial modulo f(x)
    poly->n = n;
    fmpz_init_set(poly->mod, q);
    fmpz_poly_init(poly->poly);
    fmpz_poly_init(poly->fmod);

    //Set polynomial f to x^n + 1
    fmpz_poly_set_coeff_si(poly->fmod, 0, 1);
    fmpz_poly_set_coeff_si(poly->fmod, n, 1);
}

void plwe_poly_clear(struct plwe_poly *poly) {
    poly->n = 0;
    fmpz_clear(poly->mod);
    fmpz_poly_clear(poly->poly);
    fmpz_poly_clear(poly->fmod);
}

void plwe_poly_mod_t(struct plwe_poly *poly, unsigned long t){
    fmpz_t coeff, q_2;
    fmpz_init(coeff);
    fmpz_init(q_2);

    fmpz_divexact_si(q_2, poly->mod, 2);
    unsigned long t_2 = t/2;

    for (int i = 0; i <= poly->n; i++){
        fmpz_poly_get_coeff_fmpz(coeff, poly->poly, i);

        //Correct number ranges in encrypted state
        if (fmpz_cmp(coeff, q_2) > 0) {
            fmpz_sub(coeff, coeff, poly->mod);
        }

        //Compute mod t
        fmpz_mod_ui(coeff, coeff, t);

        //Correct number ranges in decrypted state
        if (fmpz_cmp_ui(coeff, t_2) > 0) {
            fmpz_sub_ui(coeff, coeff, t);
        }

        fmpz_poly_set_coeff_fmpz(poly->poly, i, coeff);
    }

    fmpz_clear(coeff);
    fmpz_clear(q_2);
}

void plwe_poly_pmod(struct plwe_poly *poly){
    //FMod
    fmpz_poly_rem(poly->poly, poly->poly, poly->fmod);

    //Qmod
    fmpz_t coeff;
    fmpz_init(coeff);

    for (int i = 0; i <= poly->n; i++){
        fmpz_poly_get_coeff_fmpz(coeff, poly->poly, i);
        fmpz_mod(coeff, coeff, poly->mod);
        fmpz_poly_set_coeff_fmpz(poly->poly, i, coeff);
    }

    fmpz_clear(coeff);
}

inline __attribute__((always_inline)) void plwe_poly_set(struct plwe_poly *out, const struct plwe_poly *in){
    fmpz_poly_set(out->poly, in->poly);
    fmpz_poly_set(out->fmod, in->fmod);
    fmpz_set(out->mod, in->mod);
    out->n = in->n;
}

inline __attribute__((always_inline)) void plwe_poly_add(struct plwe_poly *result, const struct plwe_poly *poly1, const struct plwe_poly *poly2) {
    fmpz_poly_add(result->poly, poly1->poly, poly2->poly);
}

inline __attribute__((always_inline)) void plwe_poly_mul(struct plwe_poly *result, const struct plwe_poly *poly1, const struct plwe_poly *poly2){
    fmpz_poly_mul(result->poly, poly1->poly, poly2->poly);
}

inline __attribute__((always_inline)) void plwe_poly_scalar_mul_ui(struct plwe_poly *result, const struct plwe_poly *poly, unsigned long scalar){
    fmpz_poly_scalar_mul_ui(result->poly, poly->poly, scalar);
}

inline __attribute__((always_inline)) void plwe_poly_scalar_mul_si(struct plwe_poly *result, const struct plwe_poly *poly, signed long scalar){
    fmpz_poly_scalar_mul_si(result->poly, poly->poly, scalar);
}

inline __attribute__((always_inline)) void plwe_poly_scalar_mul_mpz(struct plwe_poly *result, const struct plwe_poly *poly, mpz_t scalar){
    fmpz_poly_scalar_mul_mpz(result->poly, poly->poly, scalar);
}

void plwe_poly_print(const struct plwe_poly *poly){
    printf("---------------------------------------------------------------------\n");
    printf("n: %ld\n", poly->n);
    printf("mod: "), fmpz_print(poly->mod), printf("\n");
    printf("fmod: "), fmpz_poly_print(poly->fmod), printf("\n");
    printf("poly: "), fmpz_poly_print(poly->poly), printf("\n");
    printf("---------------------------------------------------------------------\n");
}

void rand_poly_uniform(struct plwe_poly *poly, const unsigned long qBits) {
    mpz_t q;
    mpz_init2(q,qBits);
    for(int i = 0; i <= poly->n; i++) {
        do {
            get_random(q, qBits);
        } while (mpz_cmp_ui(q, 0) == 0);

        fmpz_t r;
        fmpz_set_mpz(r, q);
        fmpz_poly_set_coeff_fmpz(poly->poly,i,r);
    }
    mpz_clear(q);

    plwe_poly_pmod(poly);
}

void rand_poly_gauss(struct plwe_poly *poly, const double std_dev) {
    for (int i = 0; i <= poly->n; i++) {
        signed long r = (signed long) dist_gauss_ziggurat(std_dev);
        fmpz_poly_set_coeff_si(poly->poly, i, r);
    }
    plwe_poly_pmod(poly);
}
