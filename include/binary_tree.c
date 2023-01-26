#include "binary_tree.h"

#include "util.h"

#include <flint/fmpz_poly.h>

//Static declarations
/// Algorithm 1: Generating_SHE_Parameters
/// @param[out] settings Settings to save the parameters
/// @param[in] treefunc Function for tree generation of the form func(struct node *func_node, int func_rows, const int func_m[])
/// @param[in] rows Depth of the tree
/// @param[in] m Array containing the maximum column values
/// @param[in] m_len Length of m
/// @param[in] security_level Required security level
/// @param[in] improvements_factor Maximum amount of loop iterations to find a better q
static void generate_parameters(struct settings *settings, const int m[], int m_len, struct node *root, int security_level, int improvements_factor);

/// Algorithm 2: Estimate_ResultingPoly_fromArithmeticTree
/// @param[out] node Root node
/// @param[in] n Polynomial degree n
/// @param[in] b Basis b
static void estimate_poly(struct node *node, signed long n, signed int b);

/// Algorithm 3: Compute_MultDepth_fromArithmeticTree
/// @param[in] node Root node
/// @return Arithmetic depth D
static int compute_d(struct node *node);

/// Compute log_b(v) (effectively how large degree b must be to hold v using base b)
static inline __attribute__((always_inline)) int int_log(int base, int value) {
    return (int) (log(value) / log(base));
}

/// Compute the maximum of two integers
static inline __attribute__((always_inline)) int max_int(int x, int y) {
    return (x > y) ? x : y;
}

/// Compute the maximum of two unsigned long integers
static inline __attribute__((always_inline)) unsigned long max_ulong(unsigned long x, unsigned long y) {
    return (x > y) ? x : y;
}

/// Compute the minimum of two unsigned long integers
static inline __attribute__((always_inline)) unsigned long min_ulong(unsigned long x, unsigned long y) {
    return (x < y) ? x : y;
}

static void generate_parameters(struct settings *settings, const int m[], int m_len, struct node *root, int security_level, int improvements_factor) {
    //Variables
    int current_security_level;         //Actual security level of the implementation

    signed long n = 1;                  //1, not sqrt(2) because n = 2n in each iteration, not n=n^2
    unsigned long t;                    //Message space t
    int b;                              //Base b
    int D;                              //Max ciphertext length D
    double std_deviation = 8;           //Standard deviation = 8

    int max_m = 0;

    //Find largest value in m
    for(int i=0; i < m_len; i++) {
        if(m[i] > max_m) {
            max_m = m[i];
        }
    }

    mpz_t q, two_n, one, tmp;
    mpz_init(q);
    mpz_init(two_n);
    mpz_init(tmp);
    mpz_init_set_si(one, 1);

    mpf_t res, inp;
    mpf_init(res);
    mpf_init(inp);

    do {
        n = n << 1;  // n = 2n

        //Compute b
        b = max_int(ceil(pow(max_m, 1.0/(double) n)), 2) - 1;  //Start value

        do {
            b += 1;
            estimate_poly(root, n, b);
        } while (root->degree >= n);

        //Set t to next power of 2 above root->inf_norm
        t = 1 << (1 + (int) log2((int) root->inf_norm + 1));

        //Compute D
        D = 2 + compute_d(root);  //Basic length of 2 + number of multiplications

        //Compute q > 2 * l_inf * (t * std_dev * n^1,5)^D
        mpf_set_ui(res, t);                             //res = t
        mpf_set_d(inp, std_deviation);                  //inp = std_deviation
        mpf_mul(res, res, inp);                         //res = t * std_deviation
        mpf_set_d(inp, pow((double) n, 1.5 ));       //inp = n^1.5
        mpf_mul(res, res, inp);                         //res = t * std_deviation * n^1.5
        mpf_pow_ui(res, res, D + 2);                    //res = (t * std_deviation * n^1.5)^(D+2)
        mpf_mul_ui(res, res, root->inf_norm);           //res = l_inf * (t * std_deviation * n^1.5)^(D+2)
        mpf_mul_ui(res, res, 2);                        //res = 2 * l_inf * (t * std_deviation * n^1.5)^(D+2)
        mpz_set_f(q, res);                              //q = 2 * l_inf * (t * std_deviation * n^1.5)^(D+2)

        mpz_nextprime(q, q);                            //q (prime) > 2 * l_inf * (t * std_deviation * n^1.5)^(D+2)

        //Additional improvements for q
        mpz_set_si(two_n, n);                           //two_n = n
        mpz_mul_2exp(two_n, two_n, 1);                  //two_n = 2n

        int counter = 0;
        mpz_set(tmp,q);

        while ((counter < improvements_factor || improvements_factor == -1) && mpz_congruent_p(one, tmp, two_n) == 0) {
            mpz_nextprime(tmp, tmp);
            counter += 1;
        }

        if (mpz_congruent_p(one, tmp, two_n) != 0) {
            //q fulfills the requirement
            mpz_set(q, tmp);
        }

        //Compute security level = (1.8 * (2n + l)^2)/(n * qBits) - 140
        mpf_set_z(res, two_n);                                                //res = 2n
        mpf_set_si(inp, n_sizeinbase(root->inf_norm, 2));               //inp = l
        mpf_add(res, res, inp);                                               //res = 2n+l
        mpf_mul(res,res,res);                                                 //res = (2n + l)^2
        mpf_set_d(inp, 1.8);                                                  //inp = 1.8
        mpf_mul(res, res, inp);                                               //res = 1.8 * (2n + l)^2
        mpz_set_si(tmp, mpz_sizeinbase(q, 2));                                //tmp = qBits
        mpz_mul_si(tmp, tmp, n);                                              //tmp = n * qBits
        mpf_set_z(inp, tmp);                                                  //inp = n * qBits
        mpf_div(res, res, inp);                                               //res = 1.8 * (2n + l)^2 / (n * qBits)
        mpf_sub_ui(res, res, 140);                                            //res = 1.8 * (2n + l)^2 / (n * qBits) - 140

        current_security_level = mpf_get_si(res);
    } while (current_security_level < security_level);

    fmpz_t qout;
    fmpz_init(qout);
    fmpz_set_mpz(qout, q);

    settings_init(settings, n_sizeinbase(n, 2), qout, t, b, D);

    fmpz_clear(qout);

    mpf_clear(res);
    mpf_clear(inp);

    mpz_clear(one);
    mpz_clear(tmp);
    mpz_clear(two_n);
}

static void estimate_poly(struct node *node ,signed long n, signed int b) {
    if(node == NULL) {
        // Parent node is a leaf
        return;
    }

    estimate_poly(node->left_node, n, b);
    estimate_poly(node->right_node, n, b);

    if ((node->left_node == NULL) && (node->right_node == NULL)) {
        // This node is a leaf
        node->degree = int_log(b, node->M);

        if (node->M < b) {
            node->inf_norm = node->M;
        }
        else {
            node->inf_norm = b - 1;
        }
    }
    else if (node->type == plus) {
        // This node is an addition node
        node->degree = max_ulong(node->left_node->degree, node->right_node->degree);
        node->inf_norm = node->left_node->inf_norm + node->right_node->inf_norm;
    }
    else if (node->type == multiply) {
        // This node is a multiplication node
        node->degree = node->left_node->degree + node->right_node->degree;
        node->inf_norm = (2 * (min_ulong(node->left_node->degree, node->right_node->degree) + 1)) *
                         (node->left_node->inf_norm * node->right_node->inf_norm);
    }
}

static int compute_d(struct node *node) {
    if(node->type == plus) {
        return max_int(compute_d(node->left_node), compute_d(node->right_node));
    }

    if(node->type == multiply) {
        return 1 + compute_d(node->left_node) + compute_d(node->right_node);
    }

    return 0;  // return 0, it's a leaf
}

void create_tree_and_generate_params(struct settings *settings, void (*treefunc)(struct node *func_node, int func_rows, const int func_m[]), int rows, const int m[], int m_len,  int security_level, int improvements_factor) {
    //Init root
    struct node root;
    root.inf_norm = 0;
    root.degree = 0;
    root.type = plus;

    treefunc(&root, rows, m);
    generate_parameters(settings, m, m_len, &root, security_level, improvements_factor);
}
