#ifndef HEML_TRAINSVM_H_
#define HEML_TRAINSVM_H_
#include <Scheme.h>
#include <SecretKey.h>
#include <Ciphertext.h>
#include "CipherSVM.h"
#include "TimeUtils.h"

#include <complex>

using namespace std;
using namespace NTL;

class TrainSVM{
    public:
        TrainSVM(long dims, long numIters);
        
        long dim, numIter,slots, batch,sBits, pBits,wBits,bBits;
        //dim is related to the number of data in this case.
        long logQ, aBits,logN;
        Context* context;
        SecretKey* secretKey;
        Scheme* scheme;
        CipherSVM* cipherSVM;
        TimeUtils timeutils;
        //resulting encrypted vector
        Ciphertext encValw;
        Ciphertext encWData;
        //resulting decrypted vector
        double* cwtData; //w value
        double* cwtVal; //f(x) value
        
        long suggestLogN(long lambda, long logQ);
        void trainEncLGD(double** zDataTrain, double lr);
        void decAData(double* AData, Ciphertext encAData);
        void printDecCiphtxt(Ciphertext encData);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif
