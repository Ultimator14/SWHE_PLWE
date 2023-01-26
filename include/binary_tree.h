#ifndef CUSTOM_BINARY_TREE_H
#define CUSTOM_BINARY_TREE_H

//Forward declarations
struct settings;    /// defined in util.h

enum node_type {
    plus = 1,
    multiply = 2,
    value = 3,
};

struct node {
    unsigned long inf_norm;
    unsigned long degree;
    enum node_type type; // + or * or value
    int M;

    struct node *left_node;
    struct node *right_node;
};

/// Create binary tree and compute parameters
/// @param[out] settings Settings to save the parameters
/// @param[in] treefunc Function for tree generation of the form func(struct node *func_node, int func_rows, const int func_m[])
/// @param[in] rows Depth of the tree
/// @param[in] m Array containing the maximum column values
/// @param[in] m_len Length of m
/// @param[in] security_level Required security level
/// @param[in] improvements_factor Maximum amount of loop iterations to find a better q
void create_tree_and_generate_params(struct settings *settings, void (*treefunc)(struct node *func_node, int func_rows, const int func_m[]), int rows, const int m[], int m_len,  int security_level, int improvements_factor);

#endif //CUSTOM_BINARY_TREE_H
