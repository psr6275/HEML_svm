#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "TrainSVM.h"
#include "CipherSVM.h"

using namespace std;
using namespace NTL;

/*
* Function for loading A-matrix 
*/

//완성된 kernel 행렬 가져오기
/*
* what to improve?
* we read the number of instance automatically!
* we also read the original training and test matrix!
*/
double* zDataFromFileFullA(string& path, long& factorDim, long& sampleDim, long dim) { 
	double* zeData = new double[dim*dim];
/*
 dim = 64
 sampleDim = 샘플 개수 (row 개수)
 factorDim = attribute 개수 (column 개수)
 */
// test
	sampleDim = 0;	
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line, temp;
		long i;
		size_t start, end;
		while(getline(openFile, line)){
	start=0;
	long j=0;

			do {
				end = line.find_first_of (',', start);
				temp = line.substr(start,end);

				zeData[sampleDim*dim+j]=atof(temp.c_str());
				j++;
				start = end + 1;

			} while(start);

			sampleDim++;

		}
	} else {
		cout << "Error: cannot read file" << endl;
	}



	return zeData;
}



/*
* This file load the a-matrix first
* and then execute the training procedure
*/



int main(int argc, char **argv){
    
    SetNumThreads(8);
    const long dim=64;
    long factorDim = dim, sampleDim = dim;
    long numIter = 3;
    long lr = 0.1;
    string trainfile(argv[1]);
    cout<<"Start to load data"<<endl;
    double* zeData = zDataFromFileFullA(trainfile, factorDim, sampleDim, dim);
    cout<<"Start to training procedure!"<<endl;
    TrainSVM::trainEncLGD(zeData, dim, numIter,lr);

    return 0;
    
}