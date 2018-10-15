#include "TrainSVM.h"

#include "Ciphertext.h"
#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>

#include "CipherSVM.h"
#include "EvaluatorUtils.h"

long TrainSVM::suggestLogN(long lambda, long logQ){
    long NBnd = ceil(logQ * (lambda +110) /3.6);
    double logNBnd = log2((double)NBnd);
    return (long)ceil(logNBnd);
}
void TrainSVM::trainEncNLGD(double** zDataTrain, long ADim, long numIter, long lr){
    
    //initialize weights (wtData)
    double* wtData = new double[ADim];
    for(long i=0;i<dim;i++){
        wtData[i] = EvaluatorUtils::randomReal(1.0);
    }
    double* wData = new double[dim*dim];
    for(long j=0;j<dim*dim;j++){
        wData[j] = wtData[j%dim];
    }//horizontal copy generation

    //Set parameters and create scheme
    long wBits = 30;
    long pBits = 20;
    long lBits = 5;
    long aBits = 3;
    //not determine the parameters yet
    
    long sampleDim = dim, factorDim = dim;
    long kdeg = 3;
    long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDim));
    long kBits = (long)ceil(log2(kdeg));

    long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);
    long logN = TrainSVM::suggestLogN(80, logQ);
	long bBits = fdimBits;
	long batch = 1 << bBits;    
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
    cout << dim << "batch = " << batch << ", slots = " << slots  << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

    TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQ);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherSVM cipherSVM(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherSVM.generateAuxPoly(slots, batch, pBits);  //masking matrix - 1st column만 1, 나머지 0 생성
	ZZX poly2 = cipherSVM.generateAuxPoly2(slots, batch, pBits, scheme); //masking matrix - 1st row만 1, 나머지 0 생성

	timeutils.stop("Polynomial generation");



	Ciphertext encZData ;
	Ciphertext encWData ;


	timeutils.start("Encrypting Data...");
	encZData = scheme.encrypt(zDataTrain, slots, wBits, logQ);
	encWData = scheme.encrypt(wData, slots, wBits, logQ);

	timeutils.stop("Data encryption");

	timeutils.start("Precomputing");
	Ciphertext AbV= cipherSVM.GenAbVertical(encZData, scheme, poly2, bBits, wBits, pBits, slots) ; //각 column이 Ab인 행렬 생성
	Ciphertext AbH= cipherSVM.GenAbHorzon(encZData, scheme, poly, bBits, wBits, pBits, slots) ; //각 Row가 Ab인 행렬 생성
	Ciphertext AtA= cipherSVM.GenAtA(encZData, scheme, poly, poly2, bBits, wBits, pBits, batch, slots) ; //각 AtA 행렬 생성

	timeutils.stop("Precomputing Done");

	cout << " !!! START ITERATION !!! " << endl;
		
	//timeutils.start("mult iter");
    for(long iter = 0; iter < numIter; ++iter){
        cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
		cout<<"encWData.logq before: "<< encWData.logq <<endl;
        timeutils.start("Enc LGD");
        cipherSVM.encLGDiteration(AtA,AbV,AbH,encWData,poly,poly2,lr,sBits,bBits,wBits,pBits,aBits);
        timeutils.stop("Enc LGD");
        cout << "encWData.logq after: " << encWData[0].logq << endl;
        //learning 이 잘 되었는지는 어차피 Decrypt된 상태에서 하네... 일단 얘 먼저 test 해야할듯 
        }





////////// two step  : A곱하기를 두번 한 후 Ab 빼는 방법 : 걸리는 시간 약 5초
/*
 		Ciphertext encIP = encHorizonVecProduct(encZData, encWData, scheme, poly,  bBits, wBits, pBits) ;

		encWData = encVerticalVecProduct(encZData, encIP, scheme, poly2, bBits, wBits, pBits) ;
			scheme.addAndEqual(encWData, AbH);

		complex<double>* msgg = scheme.decrypt(secretKey, encWData);
*/


///////////

		////Precomputed AtA  : AtA 곱한 후 Ab빼는 방법 : 걸리는 시간 2.5초
		Ciphertext encIP = encVerticalVecProduct(AtA, encWData, scheme, poly2, bBits, wBits, pBits) ;
		scheme.addAndEqual(encIP, AbH);


		complex<double>* msgg = scheme.decrypt(secretKey, encIP);


	



		timeutils.stop("mult iter end");

	


	return 0;
}