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
        TrainSVM(long dims, long numIters);
        Scheme scheme;
        SecretKey secretKey;
        Context context;
        long dim, numIter,slots, batch, pBits;
        CipherSVM cipherSVM;
        //resulting vector
        Ciphertext encValw;
        
        long suggestLogN(long lambda, long logQ);
Scheme& 
        //id trainEncLGD(double* zDataTrain, long dim, long numIter, double lr);
        void decAlong dims, long numItersertext encAData,long wBits);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif
