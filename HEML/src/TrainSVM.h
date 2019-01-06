#ifndef HEML_TRAINSVM_H_
#define HEML_TRAINSVM_H_
#include <Scheme.h>
#include <SecretKey.h>
#include <Ciphertext.h>

#include <complex>

using namespace std;
using namespace NTL;

class TrainSVM{
    public:
        Scheme& scheme;
        SecretKey& secretKey;
        TrainSVM(Scheme& scheme, SecretKey& secretKey) : scheme(scheme), secretKey(secretKey) {}

        long suggestLogN(long lambda, long logQ);

        void trainEncLGD(double* zDataTrain, long dim, long numIter, double lr);
        void decAData(double* AData, Ciphertext encAData,long wBits);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif
