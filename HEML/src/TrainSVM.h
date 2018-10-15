#ifndef HEML_TRAINSVM_H_
#define HEML_TRAINSVM_H_

class TrainSVM{
    public:
        static long suggestLogN(long lambda, long logQ);

        static void trainEncLGD(double* zDataTrain, long dim, long numIter, double lr);
        

    
    
};
#endif