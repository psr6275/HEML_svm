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
        Scheme scheme;
        SecretKey secretKey;
        Context context;
        long dim, numIter,slots, batch, pBits,wBits,bBits;
        long logQ;
        CipherSVM cipherSVM;
        //resulting vector
        Ciphertext encValw;
        
        long suggestLogN(long lambda, long logQ);
Scheme& 
        void trainEncLGD(double* zDataTrain, double lr);
        void decAData(double* AData, Ciphertext encAData);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif
