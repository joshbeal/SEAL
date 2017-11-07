#pragma once

#include <cstdint>
#include "seal/util/defines.h"
#include "seal/util/mempool.h"
#include "seal/util/numth.h"
#include "seal/smallmodulus.h"
#include "seal/util/uintarith.h"

namespace seal
{
    namespace util
    {
        inline std::uint64_t increment_uint_mod(std::uint64_t operand, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
            if (operand >= modulus.value())
            {
                throw std::out_of_range("operand");
            }
#endif
            operand++;
            return operand - (modulus.value() & static_cast<std::uint64_t>(-static_cast<std::int64_t>(operand >= modulus.value())));
        }

        inline std::uint64_t decrement_uint_mod(std::uint64_t operand, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
            if (operand >= modulus.value())
            {
                throw std::out_of_range("operand");
            }
#endif
            std::int64_t carry = (operand == 0);
            return operand - 1 + (modulus.value() & static_cast<std::uint64_t>(-carry));
        }

        inline std::uint64_t negate_uint_mod(std::uint64_t operand, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
            if (operand >= modulus.value())
            {
                throw std::out_of_range("operand");
            }
#endif
            std::int64_t non_zero = (operand != 0);
            return (modulus.value() - operand) & static_cast<std::uint64_t>(-non_zero);
        }

        inline std::uint64_t div2_uint_mod(std::uint64_t operand, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
#endif
            if (operand & 1)
            {
                int64_t carry = add_uint64(operand, modulus.value(), 0, &operand);
                operand >>= 1;
                if (carry)
                {
                    return operand | (1ULL << (bits_per_uint64 - 1));
                }
                return operand;
            }
            return operand >> 1;
        }

        inline std::uint64_t add_uint_uint_mod(std::uint64_t operand1, std::uint64_t operand2, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
            if (operand1 >= modulus.value())
            {
                throw std::out_of_range("operand1");
            }
            if (operand2 >= modulus.value())
            {
                throw std::out_of_range("operand2");
            }
#endif
            // Sum of operands modulo SmallModulus can never wrap around 2^64
            operand1 += operand2;
            return operand1 - (modulus.value() & static_cast<std::uint64_t>(-static_cast<std::int64_t>(operand1 >= modulus.value())));
        }

        inline std::uint64_t sub_uint_uint_mod(std::uint64_t operand1, std::uint64_t operand2, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }

            if (operand1 >= modulus.value())
            {
                throw std::out_of_range("operand1");
            }
            if (operand2 >= modulus.value())
            {
                throw std::out_of_range("operand2");
            }
#endif
            std::int64_t borrow = SEAL_SUB_BORROW_UINT64(operand1, operand2, 0, &operand1);
            return operand1 + (modulus.value() & static_cast<std::uint64_t>(-borrow));
        }

        inline std::uint64_t barrett_reduce_128(const std::uint64_t *input, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (input == nullptr)
            {
                throw std::invalid_argument("input");
            }
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
#endif
            // Reduces input using base 2^64 Barrett reduction
            // input allocation size must be 128 bits

            std::uint64_t tmp1, tmp2[2], tmp3, carry;
            const std::uint64_t *const_ratio = modulus.const_ratio().data();

            // Multiply input and const_ratio
            // Round 1
            multiply_uint64_hw64(input[0], const_ratio[0], &carry);

            multiply_uint64(input[0], const_ratio[1], tmp2);
            tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, 0, &tmp1);

            // Round 2
            multiply_uint64(input[1], const_ratio[0], tmp2);
            carry = tmp2[1] + add_uint64(tmp1, tmp2[0], 0, &tmp1);

            // This is all we care about
            tmp1 = input[1] * const_ratio[1] + tmp3 + carry;

            // Barrett subtraction
            tmp3 = input[0] - tmp1 * modulus.value();

            // Claim: One more subtraction is enough
            return tmp3 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3 >= modulus.value())));
        }

        inline std::uint64_t multiply_uint_uint_mod(std::uint64_t operand1, std::uint64_t operand2, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
#endif
            std::uint64_t z[2];
            multiply_uint64(operand1, operand2, z);
            return barrett_reduce_128(z, modulus);
        }

        inline void multiply_uint_uint_mod_vector(const std::uint64_t *operand1, const std::uint64_t *operand2, int count, const SmallModulus &modulus, std::uint64_t *result)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
#endif
            int batch_count = count >> 3;
            for (int i = 0; i < batch_count << 3; i += 8)
            {
                std::uint64_t z_0[2], z_1[2], z_2[2], z_3[2],
                    z_4[2], z_5[2], z_6[2], z_7[2];

                multiply_uint64(operand1[i], operand2[i], z_0);
                multiply_uint64(operand1[i + 1], operand2[i + 1], z_1);
                multiply_uint64(operand1[i + 2], operand2[i + 2], z_2);
                multiply_uint64(operand1[i + 3], operand2[i + 3], z_3);
                multiply_uint64(operand1[i + 4], operand2[i + 4], z_4);
                multiply_uint64(operand1[i + 5], operand2[i + 5], z_5);
                multiply_uint64(operand1[i + 6], operand2[i + 6], z_6);
                multiply_uint64(operand1[i + 7], operand2[i + 7], z_7);

                const std::uint64_t *const_ratio = modulus.const_ratio().data();
                std::uint64_t carry_0, carry_1, carry_2, carry_3,
                    carry_4, carry_5, carry_6, carry_7;

                multiply_uint64_hw64(z_0[0], const_ratio[0], &carry_0);
                multiply_uint64_hw64(z_1[0], const_ratio[0], &carry_1);
                multiply_uint64_hw64(z_2[0], const_ratio[0], &carry_2);
                multiply_uint64_hw64(z_3[0], const_ratio[0], &carry_3);
                multiply_uint64_hw64(z_4[0], const_ratio[0], &carry_4);
                multiply_uint64_hw64(z_5[0], const_ratio[0], &carry_5);
                multiply_uint64_hw64(z_6[0], const_ratio[0], &carry_6);
                multiply_uint64_hw64(z_7[0], const_ratio[0], &carry_7);

                std::uint64_t tmp2_0[2], tmp2_1[2], tmp2_2[2], tmp2_3[2],
                    tmp2_4[2], tmp2_5[2], tmp2_6[2], tmp2_7[2];

                multiply_uint64(z_0[0], const_ratio[1], tmp2_0);
                multiply_uint64(z_1[0], const_ratio[1], tmp2_1);
                multiply_uint64(z_2[0], const_ratio[1], tmp2_2);
                multiply_uint64(z_3[0], const_ratio[1], tmp2_3);
                multiply_uint64(z_4[0], const_ratio[1], tmp2_4);
                multiply_uint64(z_5[0], const_ratio[1], tmp2_5);
                multiply_uint64(z_6[0], const_ratio[1], tmp2_6);
                multiply_uint64(z_7[0], const_ratio[1], tmp2_7);

                std::uint64_t tmp1_0, tmp1_1, tmp1_2, tmp1_3,
                    tmp1_4, tmp1_5, tmp1_6, tmp1_7;

                std::uint64_t tmp3_0, tmp3_1, tmp3_2, tmp3_3,
                    tmp3_4, tmp3_5, tmp3_6, tmp3_7;

                tmp3_0 = tmp2_0[1] + add_uint64(tmp2_0[0], carry_0, 0, &tmp1_0);
                tmp3_1 = tmp2_1[1] + add_uint64(tmp2_1[0], carry_1, 0, &tmp1_1);
                tmp3_2 = tmp2_2[1] + add_uint64(tmp2_2[0], carry_2, 0, &tmp1_2);
                tmp3_3 = tmp2_3[1] + add_uint64(tmp2_3[0], carry_3, 0, &tmp1_3);
                tmp3_4 = tmp2_4[1] + add_uint64(tmp2_4[0], carry_4, 0, &tmp1_4);
                tmp3_5 = tmp2_5[1] + add_uint64(tmp2_5[0], carry_5, 0, &tmp1_5);
                tmp3_6 = tmp2_6[1] + add_uint64(tmp2_6[0], carry_6, 0, &tmp1_6);
                tmp3_7 = tmp2_7[1] + add_uint64(tmp2_7[0], carry_7, 0, &tmp1_7);

                multiply_uint64(z_0[1], const_ratio[0], tmp2_0);
                multiply_uint64(z_1[1], const_ratio[0], tmp2_1);
                multiply_uint64(z_2[1], const_ratio[0], tmp2_2);
                multiply_uint64(z_3[1], const_ratio[0], tmp2_3);
                multiply_uint64(z_4[1], const_ratio[0], tmp2_4);
                multiply_uint64(z_5[1], const_ratio[0], tmp2_5);
                multiply_uint64(z_6[1], const_ratio[0], tmp2_6);
                multiply_uint64(z_7[1], const_ratio[0], tmp2_7);

                carry_0 = tmp2_0[1] + add_uint64(tmp1_0, tmp2_0[0], 0, &tmp1_0);
                carry_1 = tmp2_1[1] + add_uint64(tmp1_1, tmp2_1[0], 0, &tmp1_1);
                carry_2 = tmp2_2[1] + add_uint64(tmp1_2, tmp2_2[0], 0, &tmp1_2);
                carry_3 = tmp2_3[1] + add_uint64(tmp1_3, tmp2_3[0], 0, &tmp1_3);
                carry_4 = tmp2_4[1] + add_uint64(tmp1_4, tmp2_4[0], 0, &tmp1_4);
                carry_5 = tmp2_5[1] + add_uint64(tmp1_5, tmp2_5[0], 0, &tmp1_5);
                carry_6 = tmp2_6[1] + add_uint64(tmp1_6, tmp2_6[0], 0, &tmp1_6);
                carry_7 = tmp2_7[1] + add_uint64(tmp1_7, tmp2_7[0], 0, &tmp1_7);

                tmp1_0 = z_0[1] * const_ratio[1] + tmp3_0 + carry_0;
                tmp1_1 = z_1[1] * const_ratio[1] + tmp3_1 + carry_1;
                tmp1_2 = z_2[1] * const_ratio[1] + tmp3_2 + carry_2;
                tmp1_3 = z_3[1] * const_ratio[1] + tmp3_3 + carry_3;
                tmp1_4 = z_4[1] * const_ratio[1] + tmp3_4 + carry_4;
                tmp1_5 = z_5[1] * const_ratio[1] + tmp3_5 + carry_5;
                tmp1_6 = z_6[1] * const_ratio[1] + tmp3_6 + carry_6;
                tmp1_7 = z_7[1] * const_ratio[1] + tmp3_7 + carry_7;

                tmp3_0 = z_0[0] - tmp1_0 * modulus.value();
                tmp3_1 = z_1[0] - tmp1_1 * modulus.value();
                tmp3_2 = z_2[0] - tmp1_2 * modulus.value();
                tmp3_3 = z_3[0] - tmp1_3 * modulus.value();
                tmp3_4 = z_4[0] - tmp1_4 * modulus.value();
                tmp3_5 = z_5[0] - tmp1_5 * modulus.value();
                tmp3_6 = z_6[0] - tmp1_6 * modulus.value();
                tmp3_7 = z_7[0] - tmp1_7 * modulus.value();

                result[i] = tmp3_0 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_0 >= modulus.value())));
                result[i + 1] = tmp3_1 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_1 >= modulus.value())));
                result[i + 2] = tmp3_2 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_2 >= modulus.value())));
                result[i + 3] = tmp3_3 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_3 >= modulus.value())));
                result[i + 4] = tmp3_4 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_4 >= modulus.value())));
                result[i + 5] = tmp3_5 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_5 >= modulus.value())));
                result[i + 6] = tmp3_6 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_6 >= modulus.value())));
                result[i + 7] = tmp3_7 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_7 >= modulus.value())));
            }
            for (int i = batch_count << 3; i < count; i++)
            {
                uint64_t z[2];
                multiply_uint64(operand1[i], operand2[i], z);
                result[i] = barrett_reduce_128(z, modulus);
            }
        }

        inline void multiply_uint_scalar_mod_vector(const std::uint64_t *operand, int count, std::uint64_t scalar, const SmallModulus &modulus, std::uint64_t *result)
        {
#ifdef SEAL_DEBUG
            if (modulus.value() == 0)
            {
                throw std::invalid_argument("modulus");
            }
#endif
            int batch_count = count >> 3;
            for (int i = 0; i < batch_count << 3; i += 8)
            {
                std::uint64_t z_0[2], z_1[2], z_2[2], z_3[2],
                    z_4[2], z_5[2], z_6[2], z_7[2];

                multiply_uint64(operand[i], scalar, z_0);
                multiply_uint64(operand[i + 1], scalar, z_1);
                multiply_uint64(operand[i + 2], scalar, z_2);
                multiply_uint64(operand[i + 3], scalar, z_3);
                multiply_uint64(operand[i + 4], scalar, z_4);
                multiply_uint64(operand[i + 5], scalar, z_5);
                multiply_uint64(operand[i + 6], scalar, z_6);
                multiply_uint64(operand[i + 7], scalar, z_7);

                const std::uint64_t *const_ratio = modulus.const_ratio().data();
                std::uint64_t carry_0, carry_1, carry_2, carry_3,
                    carry_4, carry_5, carry_6, carry_7;

                multiply_uint64_hw64(z_0[0], const_ratio[0], &carry_0);
                multiply_uint64_hw64(z_1[0], const_ratio[0], &carry_1);
                multiply_uint64_hw64(z_2[0], const_ratio[0], &carry_2);
                multiply_uint64_hw64(z_3[0], const_ratio[0], &carry_3);
                multiply_uint64_hw64(z_4[0], const_ratio[0], &carry_4);
                multiply_uint64_hw64(z_5[0], const_ratio[0], &carry_5);
                multiply_uint64_hw64(z_6[0], const_ratio[0], &carry_6);
                multiply_uint64_hw64(z_7[0], const_ratio[0], &carry_7);

                std::uint64_t tmp2_0[2], tmp2_1[2], tmp2_2[2], tmp2_3[2],
                    tmp2_4[2], tmp2_5[2], tmp2_6[2], tmp2_7[2];

                multiply_uint64(z_0[0], const_ratio[1], tmp2_0);
                multiply_uint64(z_1[0], const_ratio[1], tmp2_1);
                multiply_uint64(z_2[0], const_ratio[1], tmp2_2);
                multiply_uint64(z_3[0], const_ratio[1], tmp2_3);
                multiply_uint64(z_4[0], const_ratio[1], tmp2_4);
                multiply_uint64(z_5[0], const_ratio[1], tmp2_5);
                multiply_uint64(z_6[0], const_ratio[1], tmp2_6);
                multiply_uint64(z_7[0], const_ratio[1], tmp2_7);

                std::uint64_t tmp1_0, tmp1_1, tmp1_2, tmp1_3,
                    tmp1_4, tmp1_5, tmp1_6, tmp1_7;

                std::uint64_t tmp3_0, tmp3_1, tmp3_2, tmp3_3,
                    tmp3_4, tmp3_5, tmp3_6, tmp3_7;

                tmp3_0 = tmp2_0[1] + add_uint64(tmp2_0[0], carry_0, 0, &tmp1_0);
                tmp3_1 = tmp2_1[1] + add_uint64(tmp2_1[0], carry_1, 0, &tmp1_1);
                tmp3_2 = tmp2_2[1] + add_uint64(tmp2_2[0], carry_2, 0, &tmp1_2);
                tmp3_3 = tmp2_3[1] + add_uint64(tmp2_3[0], carry_3, 0, &tmp1_3);
                tmp3_4 = tmp2_4[1] + add_uint64(tmp2_4[0], carry_4, 0, &tmp1_4);
                tmp3_5 = tmp2_5[1] + add_uint64(tmp2_5[0], carry_5, 0, &tmp1_5);
                tmp3_6 = tmp2_6[1] + add_uint64(tmp2_6[0], carry_6, 0, &tmp1_6);
                tmp3_7 = tmp2_7[1] + add_uint64(tmp2_7[0], carry_7, 0, &tmp1_7);

                multiply_uint64(z_0[1], const_ratio[0], tmp2_0);
                multiply_uint64(z_1[1], const_ratio[0], tmp2_1);
                multiply_uint64(z_2[1], const_ratio[0], tmp2_2);
                multiply_uint64(z_3[1], const_ratio[0], tmp2_3);
                multiply_uint64(z_4[1], const_ratio[0], tmp2_4);
                multiply_uint64(z_5[1], const_ratio[0], tmp2_5);
                multiply_uint64(z_6[1], const_ratio[0], tmp2_6);
                multiply_uint64(z_7[1], const_ratio[0], tmp2_7);

                carry_0 = tmp2_0[1] + add_uint64(tmp1_0, tmp2_0[0], 0, &tmp1_0);
                carry_1 = tmp2_1[1] + add_uint64(tmp1_1, tmp2_1[0], 0, &tmp1_1);
                carry_2 = tmp2_2[1] + add_uint64(tmp1_2, tmp2_2[0], 0, &tmp1_2);
                carry_3 = tmp2_3[1] + add_uint64(tmp1_3, tmp2_3[0], 0, &tmp1_3);
                carry_4 = tmp2_4[1] + add_uint64(tmp1_4, tmp2_4[0], 0, &tmp1_4);
                carry_5 = tmp2_5[1] + add_uint64(tmp1_5, tmp2_5[0], 0, &tmp1_5);
                carry_6 = tmp2_6[1] + add_uint64(tmp1_6, tmp2_6[0], 0, &tmp1_6);
                carry_7 = tmp2_7[1] + add_uint64(tmp1_7, tmp2_7[0], 0, &tmp1_7);

                tmp1_0 = z_0[1] * const_ratio[1] + tmp3_0 + carry_0;
                tmp1_1 = z_1[1] * const_ratio[1] + tmp3_1 + carry_1;
                tmp1_2 = z_2[1] * const_ratio[1] + tmp3_2 + carry_2;
                tmp1_3 = z_3[1] * const_ratio[1] + tmp3_3 + carry_3;
                tmp1_4 = z_4[1] * const_ratio[1] + tmp3_4 + carry_4;
                tmp1_5 = z_5[1] * const_ratio[1] + tmp3_5 + carry_5;
                tmp1_6 = z_6[1] * const_ratio[1] + tmp3_6 + carry_6;
                tmp1_7 = z_7[1] * const_ratio[1] + tmp3_7 + carry_7;

                tmp3_0 = z_0[0] - tmp1_0 * modulus.value();
                tmp3_1 = z_1[0] - tmp1_1 * modulus.value();
                tmp3_2 = z_2[0] - tmp1_2 * modulus.value();
                tmp3_3 = z_3[0] - tmp1_3 * modulus.value();
                tmp3_4 = z_4[0] - tmp1_4 * modulus.value();
                tmp3_5 = z_5[0] - tmp1_5 * modulus.value();
                tmp3_6 = z_6[0] - tmp1_6 * modulus.value();
                tmp3_7 = z_7[0] - tmp1_7 * modulus.value();

                result[i] = tmp3_0 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_0 >= modulus.value())));
                result[i + 1] = tmp3_1 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_1 >= modulus.value())));
                result[i + 2] = tmp3_2 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_2 >= modulus.value())));
                result[i + 3] = tmp3_3 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_3 >= modulus.value())));
                result[i + 4] = tmp3_4 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_4 >= modulus.value())));
                result[i + 5] = tmp3_5 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_5 >= modulus.value())));
                result[i + 6] = tmp3_6 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_6 >= modulus.value())));
                result[i + 7] = tmp3_7 - (modulus.value() & static_cast<uint64_t>(-static_cast<std::int64_t>(tmp3_7 >= modulus.value())));
            }
            for (int i = batch_count << 3; i < count; i++)
            {
                uint64_t z[2];
                multiply_uint64(operand[i], scalar, z);
                result[i] = barrett_reduce_128(z, modulus);
            }
        }

        inline void modulo_uint_inplace(std::uint64_t *value, int value_uint64_count, const SmallModulus &modulus)
        {
#ifdef SEAL_DEBUG
            if (value == nullptr && value_uint64_count > 0)
            {
                throw std::invalid_argument("value");
            }
            if (value_uint64_count < 0)
            {
                throw std::invalid_argument("value_uint64_count");
            }
#endif
            if (value_uint64_count == 1)
            {
                value[0] %= modulus.value();
                return;
            }

            // Starting from the top, reduce always 128-bit blocks
            for (int i = value_uint64_count - 2; i >= 0; i--)
            {
                value[i] = barrett_reduce_128(value + i, modulus);
                value[i + 1] = 0;
            }
        }

        inline std::uint64_t modulo_uint(const std::uint64_t *value, int value_uint64_count, const SmallModulus &modulus, MemoryPool &pool)
        {
#ifdef SEAL_DEBUG
            if (value == nullptr && value_uint64_count > 0)
            {
                throw std::invalid_argument("value");
            }
            if (value_uint64_count <= 0)
            {
                throw std::invalid_argument("value_uint64_count");
            }
#endif
            if (value_uint64_count == 1)
            {
                // If value < modulus no operation is needed
                return *value % modulus.value();
            }

            Pointer value_copy(allocate_uint(value_uint64_count, pool));
            set_uint_uint(value, value_uint64_count, value_copy.get());

            // Starting from the top, reduce always 128-bit blocks
            for (int i = value_uint64_count - 2; i >= 0; i--)
            {
                value_copy[i] = barrett_reduce_128(value_copy.get() + i, modulus);
            }

            return value_copy[0];
        }

        inline bool try_invert_uint_mod(uint64_t operand, const SmallModulus &modulus, std::uint64_t &result)
        {
            return try_mod_inverse(operand, modulus.value(), result);
        }

        bool is_primitive_root(std::uint64_t root, std::uint64_t degree, const SmallModulus &prime_modulus);

        // Try to find a primitive degree-th root of unity modulo small prime modulus, where degree must be a power of two.
        bool try_primitive_root(std::uint64_t degree, const SmallModulus &prime_modulus, std::uint64_t &destination, MemoryPool &pool);

        // Try to find the smallest (as integer) primitive degree-th root of unity modulo small prime modulus, where degree must be a power of two.
        bool try_minimal_primitive_root(std::uint64_t degree, const SmallModulus &prime_modulus, std::uint64_t &destination, MemoryPool &pool);
   
        std::uint64_t exponentiate_uint_mod(std::uint64_t operand, std::uint64_t exponent, const SmallModulus &modulus);
    }
}
