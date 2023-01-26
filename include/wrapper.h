#ifndef CUSTOM_WRAPPER_H
#define CUSTOM_WRAPPER_H

//Forward declarations
struct key;         /// defined in key.h
struct message;     /// defined in message.h
struct settings;    /// defined in util.h
struct plwe_poly;   /// defined in plwe_poly.h


/// Generate a random large prime of a certain bit-size and fill settings with parameters
/// @param[out] settings Empty settings
/// @param[in] n_power Power of n
/// @param[in] qBits Bit size of q
/// @param[in] t Message space
/// @param[in] b Encoding base
/// @param[in] D Maximum ciphertext length / maximum homomorphic depth minus 2
void settings_init_gen_prime(struct settings *settings, unsigned long n_power, unsigned long qBits, unsigned long t, signed int b, unsigned long D);

/// Generate a random large prime congruent mod 2n of a certain bit-size and fill settings with parameters
/// @param[out] settings Empty settings
/// @param[in] n_power Power of n
/// @param[in] qBits Bit size of q
/// @param[in] t Message space
/// @param[in] b Encoding base
/// @param[in] D Maximum ciphertext length / maximum homomorphic depth minus 2
void settings_init_gen_prime_congruent_mod_2n(struct settings *settings, unsigned long n_power, unsigned long qBits, unsigned long t, signed int b, unsigned long D);

/// Encode and encrypt a signed integer
/// @param[out] output Ciphertext
/// @param[in] input Plaintext
/// @param[in] settings Settings
/// @param[in] key Key
void encode_encrypt(struct message *output, signed long input, const struct settings *settings, const struct key *key);

/// Decrypt and decode a ciphertext
/// @param[in] input Ciphertext
/// @param[in] settings Settings
/// @param[in] key Key
/// @return Plaintext
signed int decrypt_decode(struct message *input, const struct settings *settings, const struct key *key);

/// Encode and add a signed integer to a ciphertext
/// @param[out] output Ciphertext result of the addition
/// @param[in] message Ciphertext
/// @param[in] plain Plaintext
/// @param[in] settings Settings
void encode_eval_add_plain(struct message *output, struct message *message, signed long plain, const struct settings *settings);

/// Encode and multiply a signed integer with a ciphertext
/// @param[out] output Ciphertext result of the multiplication
/// @param[in] message Ciphertext
/// @param[in] plain Plaintext
/// @param[in] settings Settings
void encode_eval_mul_plain(struct message *output, struct message *message, signed long plain, const struct settings *settings);

#endif //CUSTOM_WRAPPER_H
