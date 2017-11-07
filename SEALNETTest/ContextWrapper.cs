using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using System;
using Microsoft.Research.SEAL;

namespace SEALNETTest
{
    [TestClass]
    public class ContextWrapper
    {
        [TestMethod]
        public void ContextConstructorNET()
        {
            // Nothing set
            var parms = new EncryptionParameters();
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsFalse(qualifiers.ParametersSet);
                Assert.IsFalse(qualifiers.EnableFFT);
                Assert.IsFalse(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Not relatively prime coeff moduli
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 2, 30 });
            parms.SetPlainModulus(2);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsFalse(qualifiers.ParametersSet);
                Assert.IsFalse(qualifiers.EnableFFT);
                Assert.IsFalse(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Plain modulus not relatively prime to coeff moduli
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 17, 41 });
            parms.SetPlainModulus(34);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsFalse(qualifiers.ParametersSet);
                Assert.IsFalse(qualifiers.EnableFFT);
                Assert.IsFalse(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Plain modulus not smaller than product of coeff moduli
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 2 });
            parms.SetPlainModulus(3);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus[0].ToString() == "2");
                Assert.IsFalse(qualifiers.ParametersSet);
                Assert.IsFalse(qualifiers.EnableFFT);
                Assert.IsFalse(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // FFT poly but not NTT modulus
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 3 });
            parms.SetPlainModulus(2);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus[0].ToString() == "3");
                Assert.IsFalse(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsFalse(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Parameters OK; no fast plain lift
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 17, 41 });
            parms.SetPlainModulus(18);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus.ToString().Equals("2B9"));
                Assert.IsTrue(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsTrue(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Parameters OK; fast plain lift
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 17, 41 });
            parms.SetPlainModulus(16);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus.ToString().Equals("2B9"));
                Assert.IsTrue(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsTrue(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsTrue(qualifiers.EnableFastPlainLift);
            }

            // Parameters OK; no batching due to non-prime plain modulus
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 17, 41 });
            parms.SetPlainModulus(49);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus.ToString().Equals("2B9"));
                Assert.IsTrue(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsTrue(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Parameters OK; batching enabled
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 17, 41 });
            parms.SetPlainModulus(73);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus.ToString().Equals("2B9"));
                Assert.IsTrue(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsTrue(qualifiers.EnableNTT);
                Assert.IsTrue(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }

            // Parameters OK; batching and fast plain lift enabled
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 137, 193 });
            parms.SetPlainModulus(73);
            parms.SetNoiseStandardDeviation(3.19);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsTrue(context.TotalCoeffModulus.ToString().Equals("6749"));
                Assert.IsTrue(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsTrue(qualifiers.EnableNTT);
                Assert.IsTrue(qualifiers.EnableBatching);
                Assert.IsTrue(qualifiers.EnableFastPlainLift);
            }

            // Negative noise standard deviation
            parms.SetPolyModulus("1x^4 + 1");
            parms.SetCoeffModulus(new List<UInt64> { 137, 193 });
            parms.SetPlainModulus(73);
            parms.SetNoiseStandardDeviation(-0.1);
            {
                var context = new SEALContext(parms);
                var qualifiers = context.Qualifiers;
                Assert.IsFalse(qualifiers.ParametersSet);
                Assert.IsTrue(qualifiers.EnableFFT);
                Assert.IsFalse(qualifiers.EnableNTT);
                Assert.IsFalse(qualifiers.EnableBatching);
                Assert.IsFalse(qualifiers.EnableFastPlainLift);
            }
        }
    };
}