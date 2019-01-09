#ifndef HEML_TRAINSVM_H_
#define HEML_TRAINSVM_H_
#include <Scheme.h>
#include <SecretKey.h>
#include <Ciphertext.h>
#include "CipherSVM.h"

#include <complex>

using namespace std;
using namespace NTL;

class TrainSVM{
    public:
        TrainSVM(long dims, long numIters);
        
        long dim, numIter,slots, batch,sBits, pBits,wBits,bBits;
        long logQ, aBits,logN;
        Context context(long logN,long logQ);
        SecretKey secretKey(long logN);
        Scheme scheme(SecretKey secretKey,Context context);
        CipherSVM cipherSVM(Scheme scheme, SecretKey secretKey);
        //resulting vector
        Ciphertext encValw;
        
        long suggestLogN(long lambda, long logQ);
        void trainEncLGD(double* zDataTrain, double lr);
        void decAData(double* AData, Ciphertext encAData);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif
