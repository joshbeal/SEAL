using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.Research.SEAL;
using System.Collections.Generic;
using System;
using System.Linq;

namespace SEALNETTest
{
    [TestClass]
    public class PolyCRTBuilderWrapper
    {
        [TestMethod]
        public void BatchUnbatchVectorNET()
        {
            var parms = new EncryptionParameters();
            parms.SetPolyModulus("1x^64 + 1");
            parms.SetCoeffModulus(new List<SmallModulus> { DefaultParams.SmallMods60Bit(0) });
            parms.SetPlainModulus(257);

            var context = new SEALContext(parms);
            Assert.IsTrue(context.Qualifiers.EnableBatching);

            var crtbuilder = new PolyCRTBuilder(context);
            Assert.AreEqual(64, crtbuilder.SlotCount);
            var plain_vec = new List<UInt64>();
            for (int i = 0; i < crtbuilder.SlotCount; i++)
            {
                plain_vec.Add((UInt64)i);
            }

            var plain = new Plaintext();
            crtbuilder.Compose(plain_vec, plain);
            var plain_vec2 = new List<UInt64>();
            crtbuilder.Decompose(plain, plain_vec2);
            Assert.IsTrue(plain_vec.SequenceEqual(plain_vec2));

            for (int i = 0; i < crtbuilder.SlotCount; i++)
            {
                plain_vec[i] = 5;
            }
            crtbuilder.Compose(plain_vec, plain);
            Assert.IsTrue(plain.ToString().Equals("5"));
            crtbuilder.Decompose(plain, plain_vec2);
            Assert.IsTrue(plain_vec.SequenceEqual(plain_vec2));
        }

        [TestMethod]
        public void BatchUnbatchPlaintextNET()
        {
            var parms = new EncryptionParameters();
            parms.SetPolyModulus("1x^64 + 1");
            parms.SetCoeffModulus(new List<SmallModulus> { DefaultParams.SmallMods60Bit(0) });
            parms.SetPlainModulus(257);

            var context = new SEALContext(parms);
            Assert.IsTrue(context.Qualifiers.EnableBatching);

            var crtbuilder = new PolyCRTBuilder(context);
            Assert.AreEqual(64, crtbuilder.SlotCount);
            var plain = new Plaintext(crtbuilder.SlotCount);
            for (int i = 0; i < crtbuilder.SlotCount; i++)
            {
                plain[i] = (UInt64)i;
            }

            crtbuilder.Compose(plain);
            crtbuilder.Decompose(plain);
            for (int i = 0; i < crtbuilder.SlotCount; i++)
            {
                Assert.IsTrue(plain[i] == (UInt64)i);
            }

            for (int i = 0; i < crtbuilder.SlotCount; i++)
            {
                plain[i] = (UInt64)5;
            }
            crtbuilder.Compose(plain);
            Assert.IsTrue(plain.ToString().Equals("5"));
            crtbuilder.Decompose(plain);
            for (int i = 0; i < crtbuilder.SlotCount; i++)
            {
                Assert.IsTrue(plain[i] == (UInt64)5);
            }
        }
    }
}