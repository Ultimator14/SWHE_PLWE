#include "asym.h"

#include "key.h"
#include "message.h"
#include "plwe_poly.h"
#include "util.h"

#include <flint/fmpz_poly.h>
#include <stdio.h>

void keygen(struct key *key, const struct settings * const settings) {
    key_init(key, settings);

    //pk = (a0, b0 = a0 s + t e0)
    //sk <- gauss distribution X
    //a0 <- R_q
    //e0 <- gauss distribution X
    //b0 = a0 * s + t * e0

    plwe_poly_init(&key->sk, settings->q, settings->n);
    plwe_poly_init(&key->pk_a, settings->q, settings->n);
    plwe_poly_init(&key->pk_b, settings->q, settings->n);

    struct plwe_poly e0, a0s, te0;
    plwe_poly_init(&e0, settings->q, settings->n);
    plwe_poly_init(&a0s, settings->q, settings->n);
    plwe_poly_init(&te0, settings->q, settings->n);

    rand_poly_gauss(&key->sk, settings->std_dev);
    rand_poly_uniform(&key->pk_a, settings->qBits);
    rand_poly_gauss(&e0, settings->std_dev);

    plwe_poly_mul(&a0s, &key->pk_a, &key->sk);
    plwe_poly_scalar_mul_ui(&te0, &e0, settings->t);
    plwe_poly_add(&key->pk_b, &a0s, &te0);
    plwe_poly_pmod(&key->pk_b);

    plwe_poly_clear(&a0s);
    plwe_poly_clear(&e0);
    plwe_poly_clear(&te0);
}

void encrypt(struct message *message, const struct plwe_poly *m, const struct key *key) {
    if ((message->cIndex + 1) >= message->max_len){
        printf("Can't add more elements to the message. Maximum max_len (%ld) reached.\n"
               "Doing nothing.", message->max_len);
        return;
    }

    // Encryption
    // a = (a0*v + t*e'), a0=pka
    // b = (b0*v + t*e''), b0 = pkb
    // v, e' <- gauss distribution X
    // e'' <- gauss distribution X'

    //Init
    struct plwe_poly apub, bpub, poly1, poly2;
    plwe_poly_init(&apub, key->sk.mod, key->sk.n);        //a used for encryption
    plwe_poly_init(&bpub, key->sk.mod, key->sk.n);        //b used for encryption
    plwe_poly_init(&poly1, key->sk.mod, key->sk.n);       //working poly 1
    plwe_poly_init(&poly2, key->sk.mod, key->sk.n);       //working poly 2

    plwe_poly_init(&(message->c[0]), key->sk.mod, key->sk.n);
    plwe_poly_init(&(message->c[1]), key->sk.mod, key->sk.n);

    //Compute
    rand_poly_gauss(&poly1, key->settings.std_dev);                                //v <- dist

    plwe_poly_mul(&apub, &key->pk_a, &poly1);                                      //apub = a0*v
    plwe_poly_mul(&bpub, &key->pk_b, &poly1);                                      //bpub = b0*v

    rand_poly_gauss(&poly1, key->settings.std_dev);                                //e' <- dist
    rand_poly_gauss(&poly2, key->settings.greater_std_dev);                        //e'' <- greater dist
    plwe_poly_scalar_mul_ui(&poly1, &poly1, key->settings.t);                      //t*e'
    plwe_poly_scalar_mul_ui(&poly2, &poly2, key->settings.t);                      //t*e''

    plwe_poly_add(&apub, &apub, &poly1);                                           //apub = a0*v + t*e'
    plwe_poly_add(&bpub, &bpub, &poly2);                                           //bpub = b0*v + t*e''

    plwe_poly_add(&(message->c[0]), &bpub, m);                               //c0 = bpub + m
    plwe_poly_scalar_mul_si(&(message->c[1]), &apub, -1);              //c1 = -apub

    plwe_poly_pmod(&(message->c[0]));
    plwe_poly_pmod(&(message->c[1]));

    message->cIndex += 2;

    //Cleanup
    plwe_poly_clear(&apub);
    plwe_poly_clear(&bpub);
    plwe_poly_clear(&poly1);
    plwe_poly_clear(&poly2);
}

void encrypt_sym(struct message *message, const struct plwe_poly *m, const struct key *key) {
    if ((message->cIndex + 1) >= message->max_len){
        printf("Can't add more elements to the message. Maximum max_len (%ld) reached.\n"
               "Doing nothing.", message->max_len);
        return;
    }

    // Symmetrical Encryption
    // a = <- R_q
    // b = as + te
    //c0 = b + m
    //c1 = -a

    //Init
    struct plwe_poly poly1, poly2;
    plwe_poly_init(&poly1, key->sk.mod, key->sk.n);       //working poly 1
    plwe_poly_init(&poly2, key->sk.mod, key->sk.n);       //working poly 2
    plwe_poly_init(&(message->c[0]), key->sk.mod, key->sk.n);
    plwe_poly_init(&(message->c[1]), key->sk.mod, key->sk.n);

    //Compute
    rand_poly_gauss(&poly1, key->settings.std_dev);                                //e <- dist
    plwe_poly_scalar_mul_ui(&poly1, &poly1, key->settings.t);                      //t*e

    rand_poly_uniform(&poly2, key->settings.qBits);                                //a = <- R_q
    plwe_poly_scalar_mul_si(&(message->c[1]), &poly2, -1);            //c1 = -a

    plwe_poly_mul(&poly2, &key->sk, &poly2);                                      //a*s
    plwe_poly_add(&poly1, &poly1, &poly2);                                        //a*s + t*e
    plwe_poly_add(&(message->c[0]), &poly1, m);                             //c0 = a*s + t*e + m

    plwe_poly_pmod(&(message->c[0]));
    plwe_poly_pmod(&(message->c[1]));

    message->cIndex += 2;

    //Cleanup
    plwe_poly_clear(&poly1);
    plwe_poly_clear(&poly2);
}

void eval_add(struct message *result, struct message *message1, struct message *message2){
    //c(add) = c + c' = ((b + b'), -(a + a')) = (( a + a' )s + 2(e + e') + (m + m'), -(a + a'))
    //c0=b0v+te''+m c1=-a

    //Do padding for smaller polynomial
    unsigned int polynum1 = message1->cIndex;
    unsigned int polynum2 = message2->cIndex;
    unsigned int polynum_max = MAX(polynum1, polynum2);

    if (polynum1 < polynum2){
        unsigned int diff = polynum2 - polynum1;
        for (int i = 0; i < diff; i++){
            plwe_poly_init(&(message1->c[polynum_max - diff]), result->c->mod, result->c->n);
        }
    }
    else if (polynum1 > polynum2){
        unsigned int diff = polynum1 - polynum2;
        for (int i = 0; i < diff; i++){
            plwe_poly_init(&(message2->c[polynum_max - diff]), result->c->mod, result->c->n);
        }
    }
    //else nothing to do

    //Do computation
    struct plwe_poly *ptr = malloc(result->max_len * sizeof(struct plwe_poly));

    for (int i = 0; i < polynum_max; i++){
        plwe_poly_init(&ptr[i], message1->c[0].mod, message1->c[0].n);
        fmpz_poly_add(ptr[i].poly, message1->c[i].poly, message2->c[i].poly);
        plwe_poly_pmod(&ptr[i]);
    }

    free(result->c);
    result->c = ptr;
}

void eval_mul(struct message *result, const struct message *message1, const struct message *message2){
    //c(mult) = (c(mult, 0), c(mult, 1), c(mult, 2))
    //c(mult, 0) = c0c'0
    //c(mult, 1) = c0c'1 + c'0c1
    //c(mult, 2) = c1c'1

    if (message1->cIndex < 2 || message2->cIndex < 2) {
        printf("Error, polynomials do not match criteria!\n");
        return;
    }

    unsigned long len = message1->cIndex + message2->cIndex - 1;

    if (len > result->max_len){
        printf("Error, result message too small to hold result!\n");
        return;
    }

    //Do computations in new allocated memory and replace existing memory to prevent overwrites of data
    struct plwe_poly *ptr = malloc(result->max_len * sizeof(struct plwe_poly));

    for (int i = 0; i < len; i++){
        plwe_poly_init(&ptr[i], message1->c[0].mod, message1->c[0].n);  //Init plwe polys, take settings from message1
    }
    fmpz_poly_t temp;
    fmpz_poly_init(temp);

    for (int i = 0; i < message1->cIndex; i++){
        for(int j = 0; j < message2->cIndex; j++){
            fmpz_poly_mul(temp, message1->c[i].poly, message2->c[j].poly);  //Multiply ci * c'j
            fmpz_poly_add(ptr[i+j].poly,ptr[i+j].poly,temp);  //Group and add by index
        }
    }

    fmpz_poly_clear(temp);

    //Modulo
    for (int i = 0; i < len; i++){
        plwe_poly_pmod(&ptr[i]);
    }

    free(result->c);

    result->c = ptr;
    result->max_len = message1->max_len;
    result->cIndex = len;
}

void eval_add_plain(struct message *result, const struct message *message, const struct plwe_poly *plain) {
    plwe_poly_add(&result->c[0], &message->c[0], plain);
}

void eval_mul_plain(struct message *result, const struct message *message, const struct plwe_poly *plain) {
    for (int i = 0; i < message->cIndex; i++){
        plwe_poly_mul(&result->c[i], &message->c[i], plain);
    }
}

void decrypt(struct plwe_poly *m, const struct message *message, const struct key *key) {
    //Decryption works by calculating c_0 + c_1*s + c2*s^2 + c3*s^3 + ... + cl*s^l for l=cIndex
    struct plwe_poly powered_key, product;
    plwe_poly_init(&powered_key, key->settings.q, key->settings.n);
    plwe_poly_init(&product, key->settings.q, key->settings.n);

    plwe_poly_set(&powered_key, &key->sk);

    //Set c0
    plwe_poly_set(m, &(message->c[0]));

    //Set c1
    plwe_poly_mul(&product, &(message->c[1]), &(key->sk));  //Calculate c1*s
    plwe_poly_add(m, m, &product);//Add to previous value

    //Set c2-cl
    for (int i = 2; i < message->cIndex; i++){
        plwe_poly_mul(&powered_key, &powered_key, &(key->sk));  //Calculate s^i
        plwe_poly_mul(&product, &message->c[i], &powered_key);  //Calculate ci*s^i
        plwe_poly_add(m, m, &product);  //Add to previous value
    }

    plwe_poly_clear(&powered_key);
    plwe_poly_clear(&product);

    plwe_poly_pmod(m);
    plwe_poly_mod_t(m, key->settings.t);
}
