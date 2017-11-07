#include "CppUnitTest.h"
#include "seal/polycrt.h"
#include "seal/context.h"
#include "seal/keygenerator.h"
#include <vector>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace seal;
using namespace seal::util;
using namespace std;

namespace SEALTest
{
    TEST_CLASS(PolyCRTBuilderTest)
    {
    public:
        TEST_METHOD(BatchUnbatchMatrix)
        {
            EncryptionParameters parms;
            parms.set_poly_modulus("1x^64 + 1");
            parms.set_coeff_modulus({ small_mods_60bit(0) });
            parms.set_plain_modulus(257);

            SEALContext context(parms);
            Assert::IsTrue(context.qualifiers().enable_batching);

            PolyCRTBuilder crtbuilder(context);
            Assert::AreEqual(64, crtbuilder.slot_count());
            vector<uint64_t> plain_vec;
            for (int i = 0; i < crtbuilder.slot_count(); i++)
            {
                plain_vec.push_back(i);
            }

            Plaintext plain;
            crtbuilder.compose(plain_vec, plain);
            vector<uint64_t> plain_vec2;
            crtbuilder.decompose(plain, plain_vec2);
            Assert::IsTrue(plain_vec == plain_vec2);

            for (int i = 0; i < crtbuilder.slot_count(); i++)
            {
                plain_vec[i] = 5;
            }
            crtbuilder.compose(plain_vec, plain);
            Assert::IsTrue(plain.to_string() == "5");
            crtbuilder.decompose(plain, plain_vec2);
            Assert::IsTrue(plain_vec == plain_vec2);
        }
    };
}