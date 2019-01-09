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
        
        long dim, numIter,slots, batch, pBits,wBits,bBits;
        long logQ, aBits,logN;
        CipherSVM& cipherSVM;
        Context context(logN,logQ);
        SecretKey secretKey(logN);
        Scheme scheme(secretKey,context);
        //resulting vector
        Ciphertext encValw;
        
        long suggestLogN(long lambda, long logQ);
Scheme& 
        void trainEncLGD(double* zDataTrain, double lr);
        void decAData(double* AData, Ciphertext encAData);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif
