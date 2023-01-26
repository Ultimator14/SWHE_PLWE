#include "threading.h"

#include "asym.h"

#include <pthread.h>
#include <stdlib.h>

void * eval_add_threaded(struct eval_thread_args *struct_arg) {
    eval_add(struct_arg->result, struct_arg->message1, struct_arg->message2);
    return NULL;
}

void * eval_mul_threaded(struct eval_thread_args *struct_arg) {
    eval_mul(struct_arg->result, struct_arg->message1, struct_arg->message2);
    return NULL;
}

struct eval_thread_args * assign_thread_args(struct message *result, struct message *message1, struct message *message2) {
    struct eval_thread_args *etargs = malloc(sizeof(struct eval_thread_args));

    etargs->result = result;
    etargs->message1 = message1;
    etargs->message2 = message2;

    return etargs;
}

inline __attribute__((always_inline)) void clear_thread_args(struct eval_thread_args *etargs) {
    free(etargs);
}
