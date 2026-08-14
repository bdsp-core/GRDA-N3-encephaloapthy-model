/* macOS/BSD shim for glibc's reentrant drand48 family.
 *
 * The C++ model uses drand48_r/srand48_r, which are glibc extensions absent on
 * Darwin. The underlying generator is fully specified, so this reproduces the
 * identical number stream rather than substituting a different RNG:
 *
 *   X(n+1) = (a * X(n) + c) mod 2^48,  a = 0x5DEECE66D, c = 0xB
 *   srand48(s) sets X = (s << 16) | 0x330E
 *   drand48() returns X(n+1) / 2^48
 *
 * Included via -include so the model's own sources are left untouched.
 */
#ifndef DRAND48_R_SHIM_H
#define DRAND48_R_SHIM_H

#if defined(__APPLE__)

#include <stdint.h>

struct drand48_data {
    uint64_t __x;   /* current 48-bit state, held in the low bits */
    uint64_t __a;   /* multiplier   */
    uint64_t __c;   /* additive constant */
    int __init;
};

static inline int srand48_r(long int __seedval, struct drand48_data *__buffer) {
    if (!__buffer) return -1;
    __buffer->__x = ((uint64_t)((uint32_t)__seedval) << 16) | 0x330EULL;
    __buffer->__a = 0x5DEECE66DULL;
    __buffer->__c = 0xBULL;
    __buffer->__init = 1;
    return 0;
}

static inline int drand48_r(struct drand48_data *__buffer, double *__result) {
    if (!__buffer || !__result) return -1;
    if (!__buffer->__init) srand48_r(0, __buffer);
    __buffer->__x = (__buffer->__a * __buffer->__x + __buffer->__c)
                    & 0xFFFFFFFFFFFFULL;           /* mod 2^48 */
    *__result = (double)__buffer->__x / 281474976710656.0;  /* / 2^48 */
    return 0;
}

#endif /* __APPLE__ */
#endif /* DRAND48_R_SHIM_H */
