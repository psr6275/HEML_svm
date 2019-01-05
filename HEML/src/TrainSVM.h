#ifndef HEML_TRAINSVM_H_
#define HEML_TRAINSVM_H_

class TrainSVM{
    public:
        static long suggestLogN(long lambda, long logQ);

        static void trainEncLGD(double* zDataTrain, long dim, long numIter, double lr);
        static void decAData(double* AData, Ciphertext encAData,long wBits);
        //static void testEncLGD(double* zDataTest, bool isFirst,);
    
    
};
#endif