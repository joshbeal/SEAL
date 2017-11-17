#pragma once

// Enable extended parameter checks in Visual Studio debug mode
#if defined(_MSC_VER) && defined(_DEBUG)
#define SEAL_DEBUG
#else
// Define SEAL_DEBUG for debug mode
#undef SEAL_DEBUG
#endif

// For security reasons one should never throw when decoding fails due
// to overflow, but in some cases this might help in diagnosing problems.
#undef SEAL_THROW_ON_DECODER_OVERFLOW

// Multiplication by a plaintext zero should not be allowed, and by
// default SEAL throws an error in this case. For performance reasons
// one might want to undefine this if appropriate checks are performed
// elsewhere.
#define SEAL_THROW_ON_MULTIPLY_PLAIN_BY_ZERO

// Compile for big-endian system (not implemented)
#undef SEAL_BIG_ENDIAN

// Bound on the bit-length of user-defined moduli
#define SEAL_USER_MODULO_BIT_BOUND 60

// Bound on the number of coefficient moduli
#define SEAL_COEFF_MOD_COUNT_BOUND 62

// Maximum value for decomposition bit count
#define SEAL_DBC_MAX 60

// Minimum value for decomposition bit count
#define SEAL_DBC_MIN 1

// Debugging help
#define SEAL_ASSERT(condition) { if(!(condition)){ std::cerr << "ASSERT FAILED: "   \
    << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; } }

// Microsoft Visual Studio 2012 or newer
#if (_MSC_VER >= 1700)

// X64
#ifdef _M_X64

// Use compiler intrinsics for better performance
#define SEAL_ENABLE_INTRIN

#ifdef SEAL_ENABLE_INTRIN
#include <intrin.h>

#pragma intrinsic(_addcarry_u64)
#define SEAL_ADD_CARRY_UINT64(operand1, operand2, carry, result) _addcarry_u64(     \
    carry,                                                                          \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))

#pragma intrinsic(_subborrow_u64)
#define SEAL_SUB_BORROW_UINT64(operand1, operand2, borrow, result) _subborrow_u64(  \
    borrow,                                                                         \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))

#pragma intrinsic(_BitScanReverse64)
#define SEAL_MSB_INDEX_UINT64(result, value) _BitScanReverse64(result, value)

#pragma intrinsic(_umul128)
#define SEAL_MULTIPLY_UINT64(operand1, operand2, result128) {                       \
    result128[0] = _umul128(                                                        \
        static_cast<unsigned long long>(operand1),                                  \
        static_cast<unsigned long long>(operand2),                                  \
        reinterpret_cast<unsigned long long*>(result128 + 1));                      \
}
#define SEAL_MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                       \
    _umul128(                                                                       \
        static_cast<unsigned long long>(operand1),                                  \
        static_cast<unsigned long long>(operand2),                                  \
        reinterpret_cast<unsigned long long*>(hw64));                               \
}
#endif

#else //_M_X64

#undef SEAL_ENABLE_INTRIN

#endif //_M_X64

#endif //_MSC_VER


// GNU GCC/G++
#if defined(__GNUC__) && (__GNUC__ < 6)
#error "SEAL requires #__GNUC__ >= 6 (currently using __GNUC__)"
#endif
#if (__GNUC__ >= 6) && defined(__cplusplus)

// Read in config.h to disable unavailable intrinsics
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// Are intrinsics enabled?
#ifdef SEAL_ENABLE_INTRIN
#include <x86intrin.h>

#ifdef SEAL_ENABLE___BUILTIN_CLZLL
// Builtin
#define SEAL_MSB_INDEX_UINT64(result, value) {                                      \
    *result = 63 - __builtin_clzll(value);                                          \
}
#endif //SEAL_ENABLE___BUILTIN_CLZLL

#ifdef SEAL_ENABLE___INT128
// Builtin
#define SEAL_MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                       \
    *hw64 = static_cast<uint64_t>((static_cast<unsigned __int128>(operand1)         \
            * static_cast<unsigned __int128>(operand2)) >> 64);                     \
}
// Builtin
#define SEAL_MULTIPLY_UINT64(operand1, operand2, result128) {                       \
    unsigned __int128 product = static_cast<unsigned __int128>(operand1) * operand2;\
    result128[0] = static_cast<uint64_t>(product);                                  \
    result128[1] = product >> 64;                                                   \
}
#endif //SEAL_ENABLE___INT128

#ifdef SEAL_ENABLE__ADDCARRY_U64
#define SEAL_ADD_CARRY_UINT64(operand1, operand2, carry, result) _addcarry_u64(     \
    carry,                                                                          \
    static_cast<unsigned long long>(operand1),                                      \
    static_cast<unsigned long long>(operand2),                                      \
    reinterpret_cast<unsigned long long*>(result))
#endif //SEAL_ENABLE__ADDCARRY_U64

#ifdef SEAL_ENABLE__SUBBORROW_U64
// Warning: Note the inverted order of operand1 and operand2
#define SEAL_SUB_BORROW_UINT64(operand1, operand2, borrow, result) _subborrow_u64(  \
    borrow,                                                                         \
    static_cast<unsigned long long>(operand2),                                      \
    static_cast<unsigned long long>(operand1),                                      \
    reinterpret_cast<unsigned long long*>(result))
#endif //SEAL_ENABLE__SUBBORROW_U64

#endif //SEAL_ENABLE_INTRIN

#endif //defined(__GNUC__ >= 6) && defined(__cplusplus)

// Use generic functions as (slower) fallback
#ifndef SEAL_ADD_CARRY_UINT64
#define SEAL_ADD_CARRY_UINT64(operand1, operand2, carry, result) add_uint64_generic(operand1, operand2, carry, result)
//#pragma message("SEAL_ADD_CARRY_UINT64 not defined. Using add_uint64_generic (see util/defines.h)")
#endif

#ifndef SEAL_SUB_BORROW_UINT64
#define SEAL_SUB_BORROW_UINT64(operand1, operand2, borrow, result) sub_uint64_generic(operand1, operand2, borrow, result)
//#pragma message("SEAL_SUB_BORROW_UINT64 not defined. Using sub_uint64_generic (see util/defines.h).")
#endif

#ifndef SEAL_MULTIPLY_UINT64
#define SEAL_MULTIPLY_UINT64(operand1, operand2, result128) {                      \
    multiply_uint64_generic(operand1, operand2, result128);                        \
}
//#pragma message("SEAL_MULTIPLY_UINT64 not defined. Using multiply_uint64_generic (see util/defines.h).")
#endif

#ifndef SEAL_MULTIPLY_UINT64_HW64
#define SEAL_MULTIPLY_UINT64_HW64(operand1, operand2, hw64) {                      \
    multiply_uint64_hw64_generic(operand1, operand2, hw64);                        \
}
//#pragma message("SEAL_MULTIPLY_UINT64 not defined. Using multiply_uint64_generic (see util/defines.h).")
#endif

#ifndef SEAL_MSB_INDEX_UINT64
#define SEAL_MSB_INDEX_UINT64(result, value) get_msb_index_generic(result, value)
//#pragma message("SEAL_MSB_INDEX_UINT64 not defined. Using get_msb_index_generic (see util/defines.h).")
#endif
