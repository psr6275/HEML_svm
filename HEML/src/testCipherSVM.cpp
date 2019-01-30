#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"

#include <cmath>
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
double** zDataFromFileFullA(string& path, long& factorDim, long& sampleDim) { 
	//double** zeData = new double[dim][dim];
        vector<vector<double>> zline;
/*
 dim = 64
 sampleDim = 샘플 개수 (row 개수)
 factorDim = attribute 개수 (column 개수)
 */
// test
        //factorDim=1;
	sampleDim = 0;	
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line, temp;
		//getline(openFile,line);
                long i;
		size_t start, end;
		while(getline(openFile, line)){
	start=0;
	//long j=0;
        vector<double> vecline;

			do {
				end = line.find_first_of (',', start);
				temp = line.substr(start,end);
                                vecline.push_back(atof(temp.c_str()));
				//zeData[sampleDim][j]=atof(temp.c_str());
				//j++;
				start = end + 1;

			} while(start);
                        zline.push_back(vecline);
			sampleDim++;

		}
	} else {
		cout << "Error: cannot read file" << endl;
	}
        double** zeData = new double*[sampleDim];
        for(long j=0;j<sampleDim;++j){
            zeData[j] = new double[sampleDim];
            for(long i=0;i<sampleDim;++i){
                zeData[j][i] = zline[j][i];
            }
        }



	return zeData;
}



/*
* This file load the a-matrix first
* and then execute the training procedure
*/



int main(int argc, char **argv){
    

    SetNumThreads(8);
    /*
    const long dim=64;
    long factorDim = dim, sampleDim = dim;
    long numIter = 1;
    double lr = 1.0;
    string trainfile(argv[1]);
    cout<<"Start to load data"<<endl;
    double* zeData = zDataFromFileFullA(trainfile, factorDim, sampleDim, dim);
    cout<<"Start to training procedure!"<<endl;
	TrainSVM trainSVM(dim,numIter);
    trainSVM.trainEncLGD(zeData, lr);
    */
    //Test and Debugging
    const long dim=3;
    long factorDim = dim,sampleDim = dim;
    long numIter = 1;
    double lr = 1.0;
    
    double** zeData = new double*[dim];
    for(int i =0;i<dim;++i){
        zeData[i] = new double[dim];
    }
	zeData[0][0] = 2.0;
	zeData[0][1] = 3.0;
	zeData[0][2] = 1.0;
	zeData[1][0] = 3.0;
	zeData[1][1] = 0.0;
	zeData[1][2] = 2.0;
	zeData[2][0] = 1.0;
	zeData[2][1] = 2.0;
	zeData[2][2] = 2.0;
    //Add generate Scheme;
    long wBits = 30;
    long pBits = 20;
    long lBits = 5, aBits = 3;
    long dimBits = (long)ceil(log2(dim));
    long logQ=400;
    long logN = 14;
    long bBits = dimBits;
    long batch = 1<<bBits;
    long sBits = dimBits+bBits;
    long slots = 1<<sBits;

    Context context(logN,logQ);
    SecretKey secretKey(logN);
    Scheme scheme(secretKey,context);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    complex<double>* IMat = new complex<double>[slots];
    for(long i=0;i<slots;i++){
        IMat[i].real(0.0);
    }
    for(long i =0;i<dim;++i){
        IMat[i*batch+i].real(1.0);
    }
    cout<<"Print IMat"<<endl;
    for (long i=0;i<slots;i++){
        cout<<IMat[i]<<", ";
    }
    cout<<endl;
    Ciphertext encIMat = scheme.encrypt(IMat,slots,wBits,logQ);
    CipherSVM cipherSVM(scheme,secretKey,encIMat);
    cipherSVM.printDecCiphtxt(encIMat);
    ZZX poly = cipherSVM.generateAuxPoly(slots,batch,pBits);
    ZZX poly2 = cipherSVM.generateAuxPoly2(slots,batch,pBits);
    Ciphertext encZData;
    cipherSVM.encMatrix(encZData,zeData,dim,slots,batch,wBits,logQ);
    cout<<"Print original zeData"<<endl;
    for (long i =0;i<dim;i++){
            for(long j=0;j<dim;j++){
            cout<<zeData[i][j]<<", ";
            }
    }
    cout<<endl;
    cout<<"print ZData"<<endl;
    cipherSVM.printDecCiphtxt(encZData);
    delete[] IMat;
    Ciphertext encAtA;
    encAtA = cipherSVM.GenEncAtA(encZData,poly,poly2,bBits,wBits,pBits,batch,slots,logQ);
    cout<<"print AtA"<<endl;
    cipherSVM.printDecCiphtxt(encAtA);
    //TrainSVM trainSVM(dim,numIter);
    //trainSVM.trainEncLGD(zeData,lr);
	/*
	for (i =0;i<dim;++i){
		for(j=0;j<i;++j){
			zeData[i*dim+j] = zeData[j*dim+i];
		}
		for(j=i;j<dim;++j){
			zeData[i*dim+j] = EvaluatorUtils::randomReal(3.0);

		}
	}
	*/
    return 0;
    
}
