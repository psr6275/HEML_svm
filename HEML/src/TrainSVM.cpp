#include "TrainSVM.h"

#include <Ciphertext.h>
#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>

#include "CipherSVM.h"
#include "EvaluatorUtils.h"

/*to do in this file
1. revise suggest LogN
2. specify the proper hyperparameter for encryption (logQ..)
3. the current version is based on the pre-computed a-matrix, 
    but we need to improve this by implementing to calculate 
    the a-matrix from the original data file.
4. expandable code for split matrix storage!
    Now, we just suppose that the matrix is not that big and
    can be stored in a ciphertext. */
// Define the constructor for TrainSVM class
// The problem is that determining parameters requires number of data and factors
TrainSVM::TrainSVM(long dims, long numIters){
	dim = dims;
	numIter = numIters;

	wBits = 30;
	pBits = 20;
	long lBits = 5;
	aBits = 3;
	long dimBits =(long)ceil(log2(dim)); 
	logQ = (dimBits + wBits + lBits) + numIter * ((3 + 1) * wBits + 2 * 3 + aBits);
    long logN = TrainSVM::suggestLogN(80, logQ);
	bBits = dimBits;
	batch = 1 << bBits;    
	long sBits = dimBits + bBits;
	slots =  1 << sBits;
    cout << dim << "batch = " << batch << ", slots = " << slots  << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	//TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	context = new Context(logN, logQ);
	secretKey = new SecretKey(logN);
    scheme = new Scheme(*secretKey, *context);
	scheme->addLeftRotKeys(*secretKey);
	scheme->addRightRotKeys(*secretKey);
	timeutils.stop("Scheme generation");
	/////////// should add the identiy matrix as the input of CipherSVM
	complex<double>* IMat = new complex<double>[slots];
	for(long i = 0;i<slots;++i){
		IMat[i].real(0.0)
	}
	for(long i=0;i<dim;++i){
		IMat[i*bBits+i].real(1.0)
	}
	delete[] IMat;
	Ciphertext encIMat = scheme->encrypt(IMat, slots,wBits,logQ);

	//initialize cipherSVM: it should be shared.
	cipherSVM = new CipherSVM(*scheme, *secretKey, encIMat);
}
long TrainSVM::suggestLogN(long lambda, long logQ){
    long NBnd = ceil(logQ * (lambda +110) /3.6);
    double logNBnd = log2((double)NBnd);
    return (long)ceil(logNBnd);
}
void TrainSVM::trainEncLGD(double* zDataTrain, double lr){
    
    //zDataTrain is A matrix but neet to be vertorized! So, this shape is dim*dim.
    //dim is the number of columns or rows of A matrix!
    //initialize weights (wtData)
    double* wtData = new double[dim];
    for(long i=0;i<dim;i++){
        //wtData[i] = EvaluatorUtils::randomReal(1.0);
		wtData[i] =0.0;
    }
    double* wData = new double[dim*dim];
    for(long j=0;j<dim*dim;j++){
        wData[j] = wtData[j%dim];
    }//horizontal copy generation
	//initial vector is represented as [v1,v2,...;v1,v2,...;v1,v2,...]

	//wtData after training
	cwtData = new double[dim];
	cwtVal = new double[dim];

    //Set parameters and create scheme
    

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherSVM->generateAuxPoly(slots, batch, pBits);  //masking matrix - 1st column만 1, 나머지 0 생성
	ZZX poly2 = cipherSVM->generateAuxPoly2(slots, batch, pBits); //masking matrix - 1st row만 1, 나머지 0 생성

	timeutils.stop("Polynomial generation");



	Ciphertext encZData ;
	Ciphertext encWData ;


	timeutils.start("Encrypting Data...");
	encZData = scheme->encrypt(zDataTrain, slots, wBits, logQ);
	encWData = scheme->encrypt(wData, slots, wBits, logQ);

	timeutils.stop("Data encryption");

	timeutils.start("Precomputing");
	Ciphertext AbV= cipherSVM->GenAbVertical(encZData, poly2, bBits, wBits, pBits, slots) ; //각 column이 Ab인 행렬 생성
	Ciphertext AbH= cipherSVM->GenAbHorzon(encZData, poly, bBits, wBits, pBits, slots) ; //각 Row가 Ab인 행렬 생성
	Ciphertext AtA= cipherSVM->GenAtA(encZData, poly, poly2, bBits, wBits, pBits, batch, slots) ; //각 AtA 행렬 생성

	timeutils.stop("Precomputing Done");

	cout << " !!! START ITERATION !!! " << endl;
		
	//timeutils.start("mult iter");
    for(long iter = 0; iter < numIter; ++iter){
        cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
		cout<<"encWData.logq before: "<< encWData.logq <<endl;
        timeutils.start("Enc LGD");
        cipherSVM->encLGDiteration(AtA,AbV,AbH,encWData,poly,poly2,lr,sBits,bBits,wBits,pBits,aBits);
        timeutils.stop("Enc LGD");
        cout << "encWData.logq after: " << encWData.logq << endl;
        //learning 이 잘 되었는지는 어차피 Decrypt된 상태에서 하네... 일단 얘 먼저 test 해야할듯 
        }
	

	
	//obtain [b,ay] for testing!
	encValw = scheme->modDownTo(encZData, encWData.logq);
	scheme->multAndEqual(encValw,encWData);
	//trained enc w data and decrypted w data
	//scheme->multByPolyAndEqual(encWData,poly2,pBits);
	TrainSVM::decAData(cwtData,encWData);
	
	//comparing and testing the trained results here!
	//decrypt the trained vector
	TrainSVM::decAData(cwtVal,encValw);
	cout<<"end training"<<endl;
	cout<<"resulting W vector"<<endl;
        for(long i=0; i <dim;++i){
        cout<<cwtData[i]<<", ";
        }
        cout<<endl;

	cout<<"resulting f(x) vector"<<endl;
        for(long i=0; i <dim;++i){
        cout<<cwtVal[i]<<", ";
        }
        cout<<endl;
}	
//void TrainSVM::test
void TrainSVM::decAData(double* AData, Ciphertext encAData){
	complex<double>* dcw = scheme->decrypt(*secretKey,encAData);
	//dcw has "slot" number component
	for (long i = 0;i<dim;++i){
		AData[i] = dcw[i].real();
	}
	
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
		//Ciphertext encIP = cipherSVM.encVerticalVecProduct(AtA, encWData, scheme, poly2, bBits, wBits, pBits) ;
		//scheme.addAndEqual(encIP, AbH);


		//complex<double>* msgg = scheme.decrypt(secretKey, encIP);


	



		//timeutils.stop("mult iter end");

	



