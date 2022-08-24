/* Copyright (c) 2016, Eidgenössische Technische Hochschule Zürich and
Forschungszentrum Jülich GmbH.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#pragma once

/*
 * preprocessor macro utilities
 */

/*
 * MARS_PP_FOREACH(macro , args...)
 *   expands macro for each entry in args...
 *
 * example:
 *
 *      #define PROTO(T) T foo(T);
 *      MARS_PP_FOREACH(PROTO, int, float, double)
 *
 *  expands to
 *
 *      int foo(int); float foo(float); double foo(double);
 *
 * example:
 *
 *      #define ALLOCATE(name) int* name = new int;
 *      #define DELETE(name) delete name;
 *      #define NAMES a, b, c
 *
 *      ALLOCATE(NAMES)
 *      DELETE(NAMES)
 *
 *  expands to
 *
 *      int* a = new int; int* b = new int; int* c = new int;
 *      delete a; delete b; delete c;
 */

// Implementation macros for MARS_PP_FOREACH:

#define MARS_PP_FOREACH_1_(M, A, ...) M(A)
#define MARS_PP_FOREACH_2_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_1_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_3_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_2_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_4_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_3_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_5_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_4_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_6_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_5_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_7_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_6_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_8_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_7_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_9_(M, A, ...) \
    M(A)                              \
    MARS_PP_FOREACH_8_(M, __VA_ARGS__)
#define MARS_PP_FOREACH_10_(M, A, ...) \
    M(A)                               \
    MARS_PP_FOREACH_9_(M, __VA_ARGS__)
#define MARS_PP_GET_11TH_ARGUMENT_(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, ...) a11

// Apply macro in first argument to each of the remaining arguments (up to 10).
// Note: if __VA_ARGS__ has size N, when it is expanded the 11th argument is the MARS_PP_FOREACH_N_ macro.
#define MARS_PP_FOREACH(M, ...)                     \
    MARS_PP_GET_11TH_ARGUMENT_(__VA_ARGS__,         \
                               MARS_PP_FOREACH_10_, \
                               MARS_PP_FOREACH_9_,  \
                               MARS_PP_FOREACH_8_,  \
                               MARS_PP_FOREACH_7_,  \
                               MARS_PP_FOREACH_6_,  \
                               MARS_PP_FOREACH_5_,  \
                               MARS_PP_FOREACH_4_,  \
                               MARS_PP_FOREACH_3_,  \
                               MARS_PP_FOREACH_2_,  \
                               MARS_PP_FOREACH_1_)  \
    (M, __VA_ARGS__)
