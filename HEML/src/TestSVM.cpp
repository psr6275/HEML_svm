#include "TestSVM.h"

#include "Ciphertext.h"
#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>

#include "CipherSVM.h"
// all operations are done in plaintext domain

//void TestSVM::decWData(cwData, encWData,wBits); 
// CipherSVM has decWData
//void TestSVM::testCipherSVM
class TestSVM{
    public: 
        //void TestSVM::testPlainSVM
        static void testAccuracy(double* zDataTest, double* evalwData,long sampleDimTrain,long sampleDimTest, double& accr);
}