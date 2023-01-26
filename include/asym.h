#ifndef CUSTOM_ASYM_H
#define CUSTOM_ASYM_H

//Forward declarations
struct key;         /// defined in key.h
struct message;     /// defined in message.h
struct settings;    /// defined in util.h
struct plwe_poly;   /// defined in plwe_poly.h

/// Generate keys
/// @param[out] key Empty Key
/// @param[in] settings Settings for generation
void keygen(struct key *key, const struct settings *settings);

/// Encrypt a plaintext asymmetrically
/// @param[out] message Empty Ciphertext
/// @param[in] m Plaintext
/// @param[in] key Key for encryption (only pk is used)
void encrypt(struct message *message, const struct plwe_poly *m, const struct key *key);

/// Encrypt a plaintext symmetrically
/// @param[out] message Empty Ciphertext
/// @param[in] m Plaintext
/// @param[in] key Key for encryption (only sk is used)
void encrypt_sym(struct message *message, const struct plwe_poly *m, const struct key *key);

/// Add two ciphertexts
/// @param[out] result Result of the operation
/// @param[in] message1 Ciphertext 1
/// @param[in] message2 Ciphertext 2
void eval_add(struct message *result, struct message *message1, struct message *message2);

/// Multiply two ciphertexts
/// @param[out] result Result of the operation
/// @param[in] message1 Ciphertext 1
/// @param[in] message2 Ciphertext 2
void eval_mul(struct message *result, const struct message *message1, const struct message *message2);

/// Add a plaintext to a ciphertext
/// @param[out] result Result of the operation
/// @param[in] message Ciphertext
/// @param[in] plain Plaintext
void eval_add_plain(struct message *result, const struct message *message, const struct plwe_poly *plain);

/// Multiply a plaintext with a ciphertext
/// @param[out] result Result of the operation
/// @param[in] message Ciphertext
/// @param[in] plain Plaintext
void eval_mul_plain(struct message *result, const struct message *message, const struct plwe_poly *plain);

/// Decrypt a ciphertext
/// @param[out] m Empty Plaintext
/// @param[in] message Ciphertext
/// @param[in] key Key for decryption (only sk is used)
void decrypt(struct plwe_poly *m, const struct message *message, const struct key *key);

#endif //CUSTOM_ASYM_H
