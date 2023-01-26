//
// Asymmetric Homomorphic Encryption over PLWE
//

#include "asym.h"
#include "binary_tree.h"
#include "dist.h"
#include "key.h"
#include "message.h"
#include "threading.h"
#include "util.h"
#include "wrapper.h"

#include <pthread.h>
#include <stdbool.h>
#include <sys/time.h>

//Helper functions
static void stopwatch() {
    static bool started = false;
    static struct timeval t1, t2;

    if (started == false) {
        //Start
        gettimeofday(&t1, NULL);
        started = true;
    }
    else{
        //Stop
        gettimeofday(&t2, NULL);

        double elapsedTime;
        elapsedTime = (double)(t2.tv_sec - t1.tv_sec) * 1000.0;
        elapsedTime += (double)(t2.tv_usec - t1.tv_usec) / 1000.0;

        printf("Time: %f ms\n", elapsedTime);

        started = false;
    }
}

//Sampling
void measure_urandom_time(){
    stopwatch();

    for (int i = 0; i< 100000; i++) {
        unsigned int data[100];
        unsigned long count = 100;
        urandom(data, count);
    }

    stopwatch();
}

void measure_box_muller_time(){
    stopwatch();

    for (int i = 0; i< 1000000; i++) {
        dist_gauss_box_muller(8);
    }

    stopwatch();
}

void measure_polar_time(){
    stopwatch();

    for (int i = 0; i< 1000000; i++) {
        dist_gauss_polar(8);
    }

    stopwatch();
}

void measure_ziggurat_time(){
    stopwatch();

    for (int i = 0; i< 1000000; i++) {
        dist_gauss_ziggurat(8);
    }

    stopwatch();
}

//Encryption
void encrypt_eval_decrypt(){
    //Settings
    struct settings settings;
    settings_init_gen_prime(&settings, 14, 500, 20000, 2, 4);
    //settings_init_gen_prime_congruent_mod_2n(&settings, 14, 500, 20000, 2, 4);

    //Check settings
    settings_check(settings);

    //Keygen
    struct key key;
    keygen(&key, &settings);

    //Encrypt
    struct message enc1, enc2;
    message_init(&enc1, &settings);
    message_init(&enc2, &settings);

    encode_encrypt(&enc1, -3, &settings, &key);   //Encode and encrypt integer -3
    encode_encrypt(&enc2, 600, &settings, &key);  //Encode and encrypt integer 600

    //Eval
    eval_mul(&enc1, &enc1, &enc2);  //Compute - 3 * 600 = -1800
    eval_add(&enc1, &enc1, &enc1);  //Compute - 1800 - 1800 = - 3600

    //Decrypt
    signed int result = decrypt_decode(&enc1, &settings, &key);
    printf("Result: %d\n", result);

    //Cleanup
    message_clear(&enc1);
    message_clear(&enc2);
}

void encrypt_eval_relin_decrypt(){
    //Settings
    struct settings settings;
    settings_init_gen_prime(&settings, 10, 110, 2000, 10, 4);

    //Keygen
    printf("Keygen...\n");
    struct key key;
    keygen(&key, &settings);

    //Encrypt
    printf("Encrypt...\n");
    struct message enc1, enc2;
    message_init(&enc1, &settings);
    message_init(&enc2, &settings);

    encode_encrypt(&enc1, 2, &settings, &key);          //Encrypt integer 2
    encode_encrypt(&enc2, 40, &settings, &key);         //Encrypt integer 40

    //Eval
    printf("Eval...\n");
    eval_mul(&enc1, &enc1, &enc2);                            //Compute 2 * 40 = 80
    eval_add(&enc1, &enc1, &enc1);                            //Compute 80 + 80 = 160

    //Create eval key
    printf("Eval Keygen...\n");
    struct key_eval key_eval;
    key_init_eval(&key_eval, &key, 2);                     //Create eval key with param T=2

    //Relinearize
    printf("Relin...\n");
    message_relinearize(&enc1, &key_eval);

    //Clear eval key
    key_clear_eval(&key_eval);

    //Decrypt
    printf("Decrypt...\n");
    signed int result = decrypt_decode(&enc1, &settings, &key);
    printf("Result: %d\n", result);

    //Cleanup
    message_clear(&enc1);
    message_clear(&enc2);
}

void encrypt_eval_plain_decrypt(){
    //Settings
    struct settings settings;
    settings_init_gen_prime(&settings, 10, 100, 200000, 2, 4);

    //Keygen
    struct key key;
    keygen(&key, &settings);

    //Encrypt
    struct message enc1, enc2;
    message_init(&enc1, &settings);
    message_init(&enc2, &settings);

    encode_encrypt(&enc1, -3, &settings, &key);  //Encode and encrypt integer -3
    encode_encrypt(&enc2, 6, &settings, &key);  //Encode and encrypt integer 600

    //Eval
    eval_mul(&enc1, &enc1, &enc2);  //Compute -3 * 600 = -1800
    eval_add(&enc1, &enc1, &enc1);  //Compute -1800 - 1800 = -3600

    //Plain eval
    encode_eval_add_plain(&enc1, &enc1, 2, &settings); //Compute -3600 + 2 = -3598
    encode_eval_mul_plain(&enc1, &enc1, 2, &settings); //Compute -3598 * 2 = -7196

    //Decrypt
    signed int result = decrypt_decode(&enc1, &settings, &key);
    printf("Result: %d\n", result);

    //Cleanup
    message_clear(&enc1);
    message_clear(&enc2);
}

static void threaded_addition() {
    //Settings
    struct settings settings;
    settings_init_gen_prime(&settings, 14, 100, 20000, 2, 4);

    //Keygen
    struct key key;
    keygen(&key, &settings);

    //Encrypt
    struct message enc1, enc2, enc3, enc4;
    message_init(&enc1, &settings);
    message_init(&enc2, &settings);
    message_init(&enc3, &settings);
    message_init(&enc4, &settings);

    encode_encrypt(&enc1, 1, &settings, &key);
    encode_encrypt(&enc2, 2, &settings, &key);
    encode_encrypt(&enc3, 3, &settings, &key);
    encode_encrypt(&enc4, 4, &settings, &key);

    //Allocate memory for argument passing
    struct eval_thread_args *args1 = assign_thread_args(&enc1, &enc1, &enc2);
    struct eval_thread_args *args2 = assign_thread_args(&enc3, &enc3, &enc4);

    //Create threads
    pthread_t thread1, thread2;
    pthread_create (&thread1, NULL, (void *) eval_add_threaded, args1);   //1+2=3
    pthread_create (&thread2, NULL, (void *) eval_add_threaded, args2);   //3+4=7
    //pthread_create (&thread1, NULL, (void *) eval_mul_threaded, args1);   //1*2=2
    //pthread_create (&thread2, NULL, (void *) eval_mul_threaded, args2);   //3*4=12

    //Wait for thread finish
    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);

    eval_mul(&enc1, &enc1, &enc3);                                             //3*7=21
    eval_add(&enc1, &enc1, &enc3);                                             //2+12=14

    //Decrypt
    signed int result = decrypt_decode(&enc1, &settings, &key);
    printf("Result: %d\n", result);

    //Cleanup
    clear_thread_args(args1);
    clear_thread_args(args2);
    message_clear(&enc1);
    message_clear(&enc2);
    message_clear(&enc3);
    message_clear(&enc4);
}

void time_measurement() {
    //Settings
    struct settings settings;

    fmpz_t q;
    fmpz_init(q);

    //generate_prime(q, 64);
    //fmpz_print(q), printf("\n");

    //Use static q
    char a[] = "16556975737954184699";  //q, 64 bit
    fmpz_set_str(q, a, 10);
    settings_init(&settings, 14, q, 10000, 2, 4);
    fmpz_clear(q);

    //Keygen
    struct key key;

    printf("Keygen: ");
    stopwatch();
    for(int i = 0; i < 1; i++) {
        keygen(&key, &settings);
    }
    stopwatch();

    //Encrypt
    struct message enc1, enc2;
    message_init(&enc1, &settings);
    message_init(&enc2, &settings);

    printf("Encrypt: ");
    stopwatch();
    for(int i = 0; i < 1; i++) {
        encode_encrypt(&enc1, 1, &settings, &key);
    }
    stopwatch();

    encode_encrypt(&enc2, 1, &settings, &key);

    //Eval
    printf("Add: ");
    stopwatch();
    for(int i = 0; i < 1; i++) {
        eval_add(&enc2, &enc1, &enc1);
    }
    stopwatch();
    printf("Mul: ");
    stopwatch();
    for(int i = 0; i < 1; i++) {
        eval_mul(&enc2, &enc1, &enc1);
    }
    stopwatch();

    //Decrypt
    printf("Decrypt: ");
    stopwatch();
    signed int result;
    for(int i = 0; i < 1; i++) {
        result = decrypt_decode(&enc1, &settings, &key);
    }
    stopwatch();
    printf("Result: %d\n", result);

    //Cleanup
    message_clear(&enc1);
    message_clear(&enc2);
}

//Misc
void key_save_load(){
    struct settings s;
    settings_init_gen_prime(&s, 4, 4, 2, 2, 2);

    struct key k,k2;
    key_init(&k, &s);
    keygen(&k, &s);

    char path[] = "outfile.dat";

    plwe_poly_print(&k.sk);
    plwe_poly_print(&k.pk_a);
    plwe_poly_print(&k.pk_b);

    key_save(&k, path);
    key_load(&k2, path);

    plwe_poly_print(&k2.sk);
    plwe_poly_print(&k2.pk_a);
    plwe_poly_print(&k2.pk_b);
}

void my_treefunc(struct node *node , int rows, const int m[]) {
    //Initialize this node
    node->inf_norm = 0;
    node->degree = 0;

    if (rows > 1) {
        node->M = 0;

        //Not the last row, allocate memory for subnodes
        node->left_node = malloc(sizeof(struct node));
        node->right_node = malloc(sizeof(struct node));
    }

    if (rows > 3) {
        //In the tree
        node->type = plus;

        node->left_node->type = plus;
        node->right_node->type = plus;

        my_treefunc(node->left_node, rows - 1, m);
        my_treefunc(node->right_node, rows - 1, m);
    }
    else if (rows == 3) {
        //Third last row, add multiplication nodes
        node->type = plus;

        node->left_node->type = multiply;
        node->right_node->type = multiply;

        my_treefunc(node->left_node, rows - 1, m);
        my_treefunc(node->right_node, rows - 1, m);
    }
    else if (rows == 2) {
        //Second last row, add value nodes (leafs)
        my_treefunc(node->left_node, rows - 1, m);
        node->left_node->M = m[0];

        my_treefunc(node->left_node, rows - 1, m);
        node->right_node->M = m[1];
    }
    else if (rows == 1) {
        //Last row, add leaf, return and let parent add the actual value
        node->type = value;
        node->left_node = NULL;
        node->right_node = NULL;
    }
}

void create_params_for_sample_tree(){
    struct settings settings;
    const int m_len = 2;
    const int m[2] = {1,2};

    create_tree_and_generate_params(&settings, my_treefunc, 4, m, m_len, 128, 20);

    settings_print(settings);
}

//Main
int main() {
    ///Sampling
    //measure_urandom_time();
    //measure_box_muller_time();
    //measure_polar_time();
    //measure_ziggurat_time();

    ///Encryption
    //encrypt_eval_decrypt();
    //encrypt_eval_relin_decrypt();
    //encrypt_eval_plain_decrypt();
    //threaded_addition();
    //time_measurement();

    ///Misc
    //key_save_load();
    //create_params_for_sample_tree();

    return 0;
}
