#include "SvmGD.h"
#include "TestGD.h"

#include "Ciphertext.h"
#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include <cmath>

#include "CipherGD.h"
#include "GD.h"

//run the svm model (training and test)
void SvmGD::testEncGD(double **zDataA, double ** zDataB, long sampleDimTrain, long sampleDimTest, bool isYfirst, long numIter){
    cout << "isYfirst = " << isYfirst << ", number of iterations = " <<numIter << "train samples = "<< sampleDimTrain << ", test samples = " <<sampleDimTest<<endl;

    long fdimBits = (long)ceil(log2(factorDim)); // data encoding space? the number of tat is important ...????
    long sdimBits = (long)Ceil(log2(sampleDimTrain));

    long wBits = 30;
    long pBits = 20;
    long lBits = 5;
    long aBits = 3;
    long kBits = (long)ceil(log2(kdeg));

    //need to determine the logQ value long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits +1)*wBits +2*pBits
    long logQ = 1000;
    long logN = TestGD::suggestLogN(80,logQ);// caculate suggest LogN for SVM
    long bBits = min(logN-1-sdimBits,fdimBits);
    long batch = 1 << bBits;
    long sBits = sdimBits + bBits;
    long slots = 1 << sBits;
    long cnum = (long)ceil((double)factorDim/batch);

    cout << "scheme generating"<<endl;
    Context context(logN,logQ);
    SecretKey secretKey(logN);
    Scheme scheme(secretKey, context);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);

    CipherGD cipherGD(scheme, secretKey);


} 
