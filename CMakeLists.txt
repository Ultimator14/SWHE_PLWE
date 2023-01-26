# Required gentoo dependencies
# dev-libs/gmp
# sci-mathematics/flint
# dev-libs/libsodium (optional)

cmake_minimum_required(VERSION 3.21)
project(Custom C)

set(CMAKE_C_STANDARD 11)

# General settings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

# Flint and GMP Library
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -L/usr/lib64/libflint.so -lflint -lgmp")

# Math Library libm.so
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm")

# Libsodium library
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lsodium")

# Threading
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

# Optimizations
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")

include_directories("include")

# Add preprocessor flag for libsodium, comment to use custom implementation
add_definitions(-DLIB_SODIUM)

add_executable(Custom main.c
        include/asym.c
        include/binary_tree.c
        include/dist.c
        include/encoding.c
        include/key.c
        include/message.c
        include/plwe_poly.c
        include/threading.c
        include/util.c
        include/wrapper.c
        )

add_executable(Custom-Debug debug.c
        include/asym.c
        include/binary_tree.c
        include/dist.c
        include/encoding.c
        include/key.c
        include/message.c
        include/plwe_poly.c
        include/threading.c
        include/util.c
        include/wrapper.c
        )
