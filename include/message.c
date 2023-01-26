#include "message.h"

#include "key.h"
#include "plwe_poly.h"
#include "util.h"

#include <flint/fmpz_poly.h>

void message_init(struct message *message, const struct settings *settings) {
    message->c = (struct plwe_poly *) malloc(settings->D * sizeof(struct plwe_poly));
    message->max_len = settings->D;
    message->cIndex = 0;
}

void message_clear(struct message *message){
    free(message->c);
    message->max_len = 0;
    message->cIndex = 0;
}

void message_relinearize(struct message *message, struct key_eval *key_eval) {
    // This function takes a message with c0,c1,c2 and transforms it to a message with c0',c1'
    if (message->cIndex != 3) {
        printf("Doing nothing. Only messages with 3 ciphertext elements are supported.\n");
        return;
    }

    //Init polys
    struct plwe_poly * c2i = malloc((key_eval->l + 1) * sizeof(struct plwe_poly));  //final polynomials used to compute c_0', c_1'
    for (unsigned long i = 0; i <= key_eval->l; i++) {
        plwe_poly_init(&c2i[i], message->c[0].mod, message->c[0].n);
    }

    //Generate c2i from c2
    mpz_t coeff, t_power, remainder;
    mpz_init(coeff);
    mpz_init(t_power);
    mpz_init(remainder);

    //Convert every coefficient of c2 to the new base; note: max possible degree is signed long
    for (signed long d = 0; d <= message->c[0].n; d++) {
        fmpz_poly_get_coeff_mpz(coeff, message->c[2].poly, d);  //get d-th coefficient of c2

        //For every decimal position using the new base; note: max possible l is limited by qBits (unsigned long)
        for (unsigned long i = key_eval->l; i > 0; i--) {
            //coeff, remainder = coeff / T^i
            mpz_set_si(t_power, key_eval->T);
            mpz_pow_ui(t_power, t_power, i);
            mpz_tdiv_qr(coeff, remainder, coeff, t_power);

            //Set corresponding coeff
            fmpz_poly_set_coeff_mpz(c2i[i].poly, d, coeff);

            //Continue next round with remainder
            mpz_set(coeff, remainder);
        }

        fmpz_poly_set_coeff_mpz(c2i[0].poly, d, coeff);  //last round with T^0=1 (do this because last loop will insert maxvalue for i and therefore break pow due to mpz_t overflow
    }

    //Clear variables
    mpz_clear(coeff);
    mpz_clear(t_power);
    mpz_clear(remainder);

    //Compute new values c0', c1' using c2i
    struct plwe_poly tmp;
    plwe_poly_init(&tmp, message->c[0].mod, message->c[0].n);

    for (unsigned long i = 0; i <= key_eval->l; i++) {
        //c0'
        plwe_poly_mul(&tmp, &key_eval->ek0[i] , &c2i[i]);
        plwe_poly_add(&message->c[0], &message->c[0], &tmp);

        //c1'
        plwe_poly_mul(&tmp, &key_eval->ek1[i] , &c2i[i]);
        plwe_poly_add(&message->c[1], &message->c[1], &tmp);
    }

    plwe_poly_pmod(&message->c[0]);
    plwe_poly_pmod(&message->c[1]);

    //Clear third element of the message
    plwe_poly_clear(&message->c[2]);
    message->cIndex -= 1;

    //Cleanup
    free(c2i);
    plwe_poly_clear(&tmp);
}
