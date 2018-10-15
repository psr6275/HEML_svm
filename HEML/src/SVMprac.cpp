#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "TrainSVM.h"

using namespace std;
using namespace NTL;

/*
* This file load the a-matrix first
* and then execute the training procedure
*/
int main(int argc, char **argv){
    
    SetNumThreads(8);
    const long dim=64;
    long numIter = 3;
    long lr = 0.1;

}