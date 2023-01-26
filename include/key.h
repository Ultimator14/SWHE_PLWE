#ifndef CUSTOM_KEY_H
#define CUSTOM_KEY_H

#include "util.h"
#include "plwe_poly.h"

struct key {
    struct settings settings;
    struct plwe_poly sk; //private key
    struct plwe_poly pk_a; //public key a
    struct plwe_poly pk_b; //public key b
};

struct key_eval {
    struct plwe_poly *ek0;
    struct plwe_poly *ek1;
    unsigned long l;  //Length of the evaluation keys
    int T;  //New encoding base
};

/// Initialize a key with settings
/// @param[out] key Empty key
/// @param[in] settings Settings for initialization
void key_init(struct key *key, const struct settings *settings);

/// Initialize an evaluation key
/// @param[out] key_eval Empty evaluation key
/// @param[in] key Key
/// @param[in] T Base parameter for the evaluation key
void key_init_eval(struct key_eval *key_eval, const struct key *key, int T);

/// Clear evaluation key (free memory)
/// @param[in,out] key_eval Evaluation key
void key_clear_eval(struct key_eval *key_eval);

/// Save a key to a file
/// @param[in] key Key
/// @param[in] path Filepath
void key_save(struct key *key, const char *path);

/// Load a key from a file
/// @param[out] key Empty Key
/// @param[in] path Filepath
void key_load(struct key *key, const char *path);

#endif //CUSTOM_KEY_H
