#ifndef CUSTOM_THREADING_H
#define CUSTOM_THREADING_H

//Forward declarations
struct message;    /// defined in message.h

struct eval_thread_args {
    struct message *result;
    struct message *message1;
    struct message *message2;
};

/// Helper function to compute threaded addition using pthread; pass this to pthread_create
/// @param[in] struct_arg Arguments for the eval_add function
void * eval_add_threaded(struct eval_thread_args *struct_arg);

/// Helper function to compute threaded multiplication using pthread; pass this to pthread_create
/// @param[in] struct_arg Arguments for the eval_mul function
void * eval_mul_threaded(struct eval_thread_args *struct_arg);

/// Helper function to allocate memory and assign parameters for evaluation arguments
/// @param[out] result Result passed to either eval_add or eval_mul
/// @param[in] message1 Message1 passed to either eval_add or eval_mul
/// @param[in] message2 Message2 passed to either eval_add or eval_mul
/// @return Pointer to eval_thread_args; pass this to pthread_create
struct eval_thread_args * assign_thread_args(struct message *result, struct message *message1, struct message *message2);

/// Helper function to free memory allocated by assign_thread_args
/// @param[in] etargs Pointer to evaluation arguments
extern void clear_thread_args(struct eval_thread_args *etargs);

#endif //CUSTOM_THREADING_H
