#include "key.h"

#include "asym.h"
#include "message.h"

#include <stdio.h>

/// Save a polynomial to a file
/// @param fp File pointer
/// @param poly Polynomial
static void write_plwe_poly(struct plwe_poly *poly, FILE *fp);

/// Read a Polynomial from a file
/// @param fp File pointer
/// @param poly Polynomial
static void read_plwe_poly(struct plwe_poly *poly, FILE *fp);

void key_init(struct key *key, const struct settings *settings) {
    key->settings = *settings;
}

void key_init_eval(struct key_eval *key_eval, const struct key *key, const int T) {
    if (T < 2 || T > 62){
        printf("Error. Base values are only supported in the range 2 <= b <=62");
        return;
    }

    //Set T and l
    key_eval->T = T;
    key_eval->l = fmpz_sizeinbase(key->settings.q, T);

    //Init variables
    mpz_t t_power;
    mpz_init(t_power);

    struct message message;  //Init in each loop iteration

    struct plwe_poly m;
    plwe_poly_init(&m, key->settings.q, key->settings.n);

    //Compute ek
    key_eval->ek0 = malloc((key_eval->l + 1) * sizeof(struct plwe_poly));
    key_eval->ek1 = malloc((key_eval->l + 1) * sizeof(struct plwe_poly));

    for(unsigned long i = 0; i <= key_eval->l; i++) {
        plwe_poly_init(&key_eval->ek0[i], key->settings.q, key->settings.n);
        plwe_poly_init(&key_eval->ek1[i], key->settings.q, key->settings.n);

        plwe_poly_mul(&m, &key->sk, &key->sk);                             //te = s^2
        mpz_set_si(t_power, T);
        mpz_pow_ui(t_power, t_power, i);
        plwe_poly_scalar_mul_mpz(&m, &m, t_power);                        //te = s^2 * T^i

        //"Encrypt" T^i * s^2
        message_init(&message, &key->settings);
        encrypt_sym(&message, &m, key);

        //Copy result to ek
        plwe_poly_set(&key_eval->ek0[i], &message.c[0]);
        plwe_poly_set(&key_eval->ek1[i], &message.c[1]);

        //Cleanup
        message_clear(&message);
    }

    mpz_clear(t_power);
}

void key_clear_eval(struct key_eval *key_eval) {
    key_eval->T = 0;
    key_eval->l = 0;
    free(key_eval->ek0);
    free(key_eval->ek1);
}

static void write_plwe_poly(struct plwe_poly *poly, FILE *fp) {
    fprintf(fp, " %ld ", poly->n);       //n
    fmpz_out_raw(fp, poly->mod);                //q
    fmpz_poly_fprint(fp, poly->poly);           //poly
}

static void read_plwe_poly(struct plwe_poly *poly, FILE *fp) {
    long n;
    fmpz_t q;
    fmpz_init(q);

    int result = fscanf(fp, " %ld ", &n);             //n

    if (result != 1) {
        //return value is number of matched elements
        printf("Error, fscanf failed!");
    }

    fmpz_inp_raw(q, fp);                                     //q

    plwe_poly_init(poly, q, n);

    fmpz_clear(q);

    fmpz_poly_fread(fp, poly->poly);                         //poly
}

void key_save(struct key *key, const char *path) {
    FILE *fp;
    fp = fopen(path, "w");

    write_plwe_poly(&key->sk, fp);
    write_plwe_poly(&key->pk_a, fp);
    write_plwe_poly(&key->pk_b, fp);

    fclose(fp);
}

void key_load(struct key *key, const char *path) {
    FILE *fp;
    fp = fopen(path, "r");

    read_plwe_poly(&key->sk, fp);
    read_plwe_poly(&key->pk_a, fp);
    read_plwe_poly(&key->pk_b, fp);

    fclose(fp);
}
