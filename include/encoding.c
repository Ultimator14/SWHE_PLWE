#include "encoding.h"

#include "plwe_poly.h"

void encode(struct plwe_poly *output, const mpz_t input, const signed int b){
    if (b < 2 || b > 62){
        printf("Error. Base values are only supported in the range 2 <= b <=62");
        return;
    }

    mpz_t scalar, remainder;
    mpz_init(scalar);
    mpz_init(remainder);

    size_t size = mpz_sizeinbase(input, b);
    mpz_set(scalar, input);

    for (int i = 0; i < size; i++){
        mpz_tdiv_qr_ui(scalar, remainder, scalar, b);
        fmpz_poly_set_coeff_mpz(output->poly, i, remainder);
    }

    mpz_clear(scalar);
    mpz_clear(remainder);
}

void decode(mpz_t output, const struct plwe_poly *input, const signed int b){
    if (b < 2 || b > 62){
        printf("Error. Base values are only supported in the range 2 <= b <=62");
        return;
    }

    mpz_t coeff, base;
    mpz_init(coeff);
    mpz_init_set_si(base, 1);

    for (int i = 0; i < input->n; i++){
        fmpz_poly_get_coeff_mpz(coeff, input->poly, i);
        mpz_mul(coeff, coeff, base);
        mpz_add(output, output, coeff);
        mpz_mul_si(base, base, b);
    }

    mpz_clear(coeff);
    mpz_clear(base);
}
