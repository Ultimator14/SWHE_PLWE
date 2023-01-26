#ifndef CUSTOM_MESSAGE_H
#define CUSTOM_MESSAGE_H

//Forward declarations
struct key_eval;    /// defined in key.h
struct settings;    /// defined in util.h

struct message {
    struct plwe_poly *c;
    unsigned long max_len;
    unsigned long cIndex;
};

/// Initialize a ciphertext
/// @param message[out] Ciphertext to initialize
/// @param settings[in] Settings for initialization
void message_init(struct message *message, const struct settings *settings);

/// Clear Ciphertext
/// @param message[in] Ciphertext
void message_clear(struct message *message);

/// Relinearize a ciphertext (reduce its elements by one)
/// @param message Ciphertext
/// @param key_eval Evaluation Key
void message_relinearize(struct message *message, struct key_eval *key_eval);

#endif //CUSTOM_MESSAGE_H
