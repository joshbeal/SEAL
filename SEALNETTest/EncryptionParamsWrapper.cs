using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.IO;
using System;
using System.Collections.Generic;
using Microsoft.Research.SEAL;
using System.Linq;

namespace SEALNETTest
{
    [TestClass]
    public class EncryptionParamsWrapper
    {
        [TestMethod]
        public void EncryptionParametersSetNET()
        {
            var parms = new EncryptionParameters();
            parms.SetNoiseStandardDeviation(0);
            parms.SetCoeffModulus(new List<SmallModulus> { });
            parms.SetPlainModulus(0);
            parms.SetPolyModulus("0");

            Assert.AreEqual(0.0, parms.NoiseStandardDeviation);
            Assert.AreEqual(0.0, parms.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus.Count == 0);
            Assert.IsTrue(parms.PlainModulus.Value == 0);
            Assert.IsTrue(parms.PolyModulus.ToString() == "0");

            parms.SetNoiseStandardDeviation(0);
            parms.SetCoeffModulus(new List<UInt64> { 2 });
            parms.SetPlainModulus(2);
            parms.SetPolyModulus("1");

            Assert.AreEqual(0.0, parms.NoiseStandardDeviation);
            Assert.AreEqual(0.0, parms.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus[0].Value == 2);
            Assert.IsTrue(parms.PlainModulus.Value == 2);
            Assert.IsTrue(parms.PolyModulus.ToString() == "1");

            parms.SetNoiseStandardDeviation(3.19);
            parms.SetCoeffModulus(new List<UInt64> { 2, 3 });
            parms.SetPlainModulus(2);
            parms.SetPolyModulus("1x^2 + 1");

            Assert.AreEqual(3.19, parms.NoiseStandardDeviation);
            Assert.AreEqual(3.19 * 6, parms.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus[0].Value == 2);
            Assert.IsTrue(parms.CoeffModulus[1].Value == 3);
            Assert.IsTrue(parms.PlainModulus.Value == 2);
            Assert.IsTrue(parms.PolyModulus.ToString() == "1x^2 + 1");

            parms.SetNoiseStandardDeviation(3.19);
            parms.SetCoeffModulus(new List<SmallModulus> { 
                DefaultParams.SmallMods30Bit(0),
                DefaultParams.SmallMods40Bit(0),
                DefaultParams.SmallMods50Bit(0) }
            );
            parms.SetPlainModulus(2);
            parms.SetPolyModulus("1x^128 + 1");

            Assert.AreEqual(3.19, parms.NoiseStandardDeviation);
            Assert.AreEqual(3.19 * 6, parms.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus[0].Equals(DefaultParams.SmallMods30Bit(0)));
            Assert.IsTrue(parms.CoeffModulus[1].Equals(DefaultParams.SmallMods40Bit(0)));
            Assert.IsTrue(parms.CoeffModulus[2].Equals(DefaultParams.SmallMods50Bit(0)));
            Assert.IsTrue(parms.PolyModulus.ToString() == "1x^128 + 1");
        }

        [TestMethod]
        public void EncryptionParametersCompareNET()
        {
            var parms1 = new EncryptionParameters();
            parms1.SetNoiseStandardDeviation(3.19);
            parms1.SetCoeffModulus(new List<SmallModulus> { DefaultParams.SmallMods30Bit(0) });
            parms1.SetPlainModulus(1 << 6);
            parms1.SetPolyModulus("1x^64 + 1");

            var parms2 = new EncryptionParameters(parms1);
            Assert.IsTrue(parms1.Equals(parms2));

            var parms3 = new EncryptionParameters();
            parms3.Set(parms2);
            Assert.IsTrue(parms3.Equals(parms2));
            parms3.SetCoeffModulus(new List<SmallModulus>{ DefaultParams.SmallMods30Bit(1) });
            Assert.IsFalse(parms3.Equals(parms2));

            parms3.Set(parms2);
            Assert.IsTrue(parms3.Equals(parms2));
            parms3.SetCoeffModulus(new List<SmallModulus> {
                DefaultParams.SmallMods30Bit(0), DefaultParams.SmallMods30Bit(1)
            });
            Assert.IsFalse(parms3.Equals(parms2));

            parms3.Set(parms2);
            parms3.SetPolyModulus("1x^128 + 1");
            Assert.IsFalse(parms3.Equals(parms1));

            parms3.Set(parms2);
            parms3.SetPlainModulus((1 << 6) + 1);
            Assert.IsFalse(parms3.Equals(parms2));

            parms3.Set(parms2);
            parms3.SetNoiseStandardDeviation(3.18);
            Assert.IsFalse(parms3.Equals(parms2));

            parms3.Set(parms2);
            parms3.SetPolyModulus("1");
            parms3.SetPolyModulus("1x^128 + 1");
            parms3.SetPolyModulus("1x^64 + 1");
            Assert.IsTrue(parms3.Equals(parms1));

            parms3.Set(parms2);
            parms3.SetCoeffModulus(new List<UInt64> { 2 });
            parms3.SetCoeffModulus(new List<SmallModulus> { DefaultParams.SmallMods50Bit(0) });
            parms3.SetCoeffModulus(parms2.CoeffModulus);
            Assert.IsTrue(parms3.Equals(parms2));
        }

        [TestMethod]
        public void EncryptionParametersSaveLoadNET()
        {
            var stream = new MemoryStream();
            var parms = new EncryptionParameters();
            var parms2 = new EncryptionParameters();
            parms.Save(stream);
            stream.Seek(0, SeekOrigin.Begin);
            parms2.Load(stream);
            Assert.AreEqual(parms.NoiseStandardDeviation, parms2.NoiseStandardDeviation);
            Assert.AreEqual(parms.NoiseMaxDeviation, parms2.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus.SequenceEqual(parms2.CoeffModulus));
            Assert.IsTrue(parms.PlainModulus.Equals(parms2.PlainModulus));
            Assert.IsTrue(parms.PolyModulus.Equals(parms2.PolyModulus));
            Assert.IsTrue(parms.Equals(parms2));

            parms.SetNoiseStandardDeviation(3.19);
            parms.SetCoeffModulus(new List<SmallModulus> { DefaultParams.SmallMods30Bit(0) });
            parms.SetPlainModulus(1 << 6);
            parms.SetPolyModulus("1x^64 + 1");
            stream.Seek(0, SeekOrigin.Begin);
            parms.Save(stream);
            stream.Seek(0, SeekOrigin.Begin);
            parms2.Load(stream);
            Assert.AreEqual(parms.NoiseStandardDeviation, parms2.NoiseStandardDeviation);
            Assert.AreEqual(parms.NoiseMaxDeviation, parms2.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus.SequenceEqual(parms2.CoeffModulus));
            Assert.IsTrue(parms.PlainModulus.Equals(parms2.PlainModulus));
            Assert.IsTrue(parms.PolyModulus.Equals(parms2.PolyModulus));
            Assert.IsTrue(parms.Equals(parms2));

            parms.SetNoiseStandardDeviation(3.19);
            parms.SetCoeffModulus(new List<SmallModulus> {
                DefaultParams.SmallMods30Bit(0),
                DefaultParams.SmallMods60Bit(0),
                DefaultParams.SmallMods60Bit(1)
            });
            parms.SetPlainModulus(1 << 30);
            parms.SetPolyModulus("1x^256 + 1");
            stream.Seek(0, SeekOrigin.Begin);
            parms.Save(stream);
            stream.Seek(0, SeekOrigin.Begin);
            parms2.Load(stream);
            Assert.AreEqual(parms.NoiseStandardDeviation, parms2.NoiseStandardDeviation);
            Assert.AreEqual(parms.NoiseMaxDeviation, parms2.NoiseMaxDeviation);
            Assert.IsTrue(parms.CoeffModulus.SequenceEqual(parms2.CoeffModulus));
            Assert.IsTrue(parms.PlainModulus.Equals(parms2.PlainModulus));
            Assert.IsTrue(parms.PolyModulus.Equals(parms2.PolyModulus));
            Assert.IsTrue(parms.Equals(parms2));
        }
    }
}