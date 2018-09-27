#include <iostream>


#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include "NTL/ZZX.h"
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <cmath>


#include "CipherGD.h"
#include "TestGD.h"

#include "GD.h"
#include "Common.h"
#include "Ciphertext.h"
#include "EvaluatorUtils.h"
#include "NumUtils.h"
#include "Scheme.h"
#include "SchemeAlgo.h"
#include "SecretKey.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "Context.h"

using namespace std;
using namespace NTL;

Ciphertext GenAtA(Ciphertext encZData, Scheme& scheme, ZZX& poly, ZZX& poly2, long bBits, long wBits, long pBits, long batch, long slots) {//AtA만드는 함수

	
	Ciphertext AtA = scheme.multByPoly(encZData, poly2, pBits);
		for (long l = 0; l < bBits; ++l) {
			Ciphertext tmp = scheme.rightRotateByPo2(AtA, l+bBits);//
			scheme.addAndEqual(AtA, tmp);
		}

		scheme.multAndEqual(AtA, encZData); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(AtA, bBits+l);//paratmeter check
			scheme.addAndEqual(AtA, rot);
		}
	
	scheme.reScaleByAndEqual(AtA, wBits);
	

	scheme.multByPolyAndEqual(AtA, poly2, pBits);
		
/////////////
		for(int i=1; i<batch; i++){
				cout << i  << "AtA row gen" << endl;

		Ciphertext encIPvec = scheme.leftRotate(encZData, (batch*i));
		scheme.multByPolyAndEqual(encIPvec, poly2, pBits);
		for (long l = 0; l < bBits; ++l) {
			Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l+bBits);//parameter check
			scheme.addAndEqual(encIPvec, tmp);
		}

		scheme.multAndEqual(encIPvec, encZData); // xy * w
		for (long l = 0; l < bBits; ++l) {

			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, bBits+l);//paratmeter check
			scheme.addAndEqual(encIPvec, rot);
		}
	
	scheme.reScaleByAndEqual(encIPvec, wBits);

				cout << i  << "AtA row gen ff" << endl;

	scheme.multByPolyAndEqual(encIPvec, poly2, pBits);
	encIPvec = scheme.rightRotate(encIPvec, batch*i);
	scheme.addAndEqual(AtA, encIPvec);
}


	return AtA;
}

//각 column이 Ab인 행렬 생성
Ciphertext GenAbHorzon(Ciphertext encZData, Scheme& scheme, ZZX& poly, long bBits, long wBits, long pBits, long slots) {
	
	Ciphertext encIPvec=encZData;
				cout << "AbH" << endl;

	complex<double>* pvals = new complex<double>[slots];
		for (long j = 0; j < slots; j++) {   //parameter check
		pvals[j].real(1.0);
	}
	ZZX ptmp = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;

	ptmp = ptmp - poly;  
	scheme.multByPolyAndEqual(encIPvec, ptmp, pBits);
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, l);
			scheme.addAndEqual(encIPvec, rot);
		}


	scheme.reScaleByAndEqual(encIPvec, wBits);
	scheme.multByPolyAndEqual(encIPvec, poly, pBits);
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l);
		scheme.addAndEqual(encIPvec, tmp);
	}

	scheme.negate(encIPvec);

	return encIPvec;
}


//각 row가 Ab인 행렬 생성
Ciphertext GenAbVertical(Ciphertext encZData, Scheme& scheme, ZZX& poly, long bBits, long wBits, long pBits, long slots) {
		
	Ciphertext encIPvec=encZData;

		complex<double>* pvals = new complex<double>[slots];
		for (long j = (1<<bBits); j < slots; j++) {   //parameter check
		pvals[j].real(1.0);
	}
	ZZX ptmp = scheme.context.encode(pvals, slots, pBits);

	delete[] pvals;
	scheme.multByPolyAndEqual(encIPvec, ptmp, pBits);

	for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, bBits+l);//paratmeter check
			scheme.addAndEqual(encIPvec, rot);
		}

	scheme.reScaleByAndEqual(encIPvec, wBits);
	
		scheme.multByPolyAndEqual(encIPvec, poly, pBits); //Vert poly
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l+bBits);//parameter check
		scheme.addAndEqual(encIPvec, tmp);
	}

		scheme.negate(encIPvec);
	return encIPvec;
}



//한 벡터의 copy가 각 row를 이루는 행렬을 입력으로 할 때, 출력값은 행렬*벡터 결과값을 한 column로 하고, 모든 column이 그 copy인 행렬인 함수

Ciphertext encHorizonVecProduct(Ciphertext encZData, Ciphertext encWData, Scheme& scheme, ZZX& poly, long bBits, long wBits, long pBits) {
	Ciphertext encIPvec;
		encIPvec = scheme.modDownTo(encZData, encWData.logq);
		scheme.multAndEqual(encIPvec, encWData); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, l);
			scheme.addAndEqual(encIPvec, rot);
		}
	
	

	scheme.reScaleByAndEqual(encIPvec, wBits);
	
	scheme.multByPolyAndEqual(encIPvec, poly, pBits);
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l);
		scheme.addAndEqual(encIPvec, tmp);
	}
	return encIPvec;
}

//한 벡터의 copy가 각 column을 이루는 행렬을 입력으로 할 때, 출력값은 행렬*벡터 결과값을 한 row로 하고, 모든 row가 그 copy인 행렬인 함수
Ciphertext encVerticalVecProduct(Ciphertext encZData, Ciphertext encWData, Scheme& scheme, ZZX& poly,  long bBits, long wBits, long pBits) {
	Ciphertext encIPvec;
		encIPvec = scheme.modDownTo(encZData, encWData.logq);
		scheme.multAndEqual(encIPvec, encWData); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, bBits+l);// 
			scheme.addAndEqual(encIPvec, rot);
		}
	
	

	scheme.reScaleByAndEqual(encIPvec, wBits);
	
	scheme.multByPolyAndEqual(encIPvec, poly, pBits); //Vert poly
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l+bBits);//parameter check
		scheme.addAndEqual(encIPvec, tmp);
	}
	return encIPvec;
}

//마스킹 행렬 - 첫 row만 1 이고 나머지 0인 행렬 만들기
ZZX generateAuxPoly2(long slots, long batch, long pBits, Scheme& scheme) {
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < batch; j ++) {   //parameter check
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}

//평문상태로 행렬-벡터 곱셈
double* rawmult(double** zData, double* wtData, long dim){
		double* vdata = new double[dim];
		for(int i=0; i< dim; i++)vdata[i]=0;

		for(int i = 0; i < dim; i++){
			for(int j = 0 ; j < dim; j++){
				vdata[j]+=zData[i][j]*wtData[j];
			}
		}

		return vdata;

}


//일반 행렬 가져오기 - GD 에 있는 함수 복붙. 우리 상황에 맞게 조정 예정
double** zDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst) {
	vector<vector < double > > zline;
	factorDim = 1; 	// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line, temp;
		getline(openFile, line);
		long i;
		size_t start, end;
		for(i = 0; i < line.length(); ++i) if(line[i] == ',' ) factorDim++;

		while(getline(openFile, line)){
			vector < double > vecline;
			do {
				end = line.find_first_of (',', start);
				temp = line.substr(start,end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while(start);
			zline.push_back(vecline);
			sampleDim++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}

	double** zData = new double*[sampleDim];
	if(isfirst) {
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			zj[0] = 2 * zline[j][0] - 1;
			for(long i = 1; i < factorDim; ++i){
				zj[i] = zj[0] * zline[j][i];
			}
			zData[j] = zj;
		}
	} else {
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			zj[0] = 2 * zline[j][factorDim - 1] - 1;
			for(long i = 1; i < factorDim; ++i){
				zj[i] = zj[0] * zline[j][i-1];
			}
			zData[j] = zj;
		}
	}
	return zData;
}



double* zDataFromFileFullA(string& path, long& factorDim, long& sampleDim, long dim) { //완성된 kernel 행렬 가져오기
	double* zeData = new double[dim*dim];

	sampleDim = 0;	// number of samples
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

//미완
Ciphertext makeMatrixA(Ciphertext enczData, Scheme& scheme, long bBits, ZZX& poly, long pBits, long wBits){
		Ciphertext Row;
		Row = scheme.rightRotate(enczData, 1);
		//첫 슬롯에 1 더하기
		for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(Row, l+bBits);// 
		scheme.addAndEqual(Row, tmp);
		}

		Ciphertext Column;
		Ciphertext tmp;
		tmp = scheme.rightRotate(enczData, 1);
		for (long l = 1; l < (1<<bBits); ++l) {
		tmp = scheme.rightRotateByPo2(tmp, bBits);
		tmp = scheme.leftRotate(tmp, 1);
		scheme.addAndEqual(Column, tmp);
	}

		scheme.multByPolyAndEqual(Column, poly, pBits); 
		for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(Column, l);
		scheme.addAndEqual(Column, tmp);
	}
	


///////////////

	return Column;


}


//Ciphertext makeMatrix(Ciphertext enczData, long& factorDim, long& sampleDim, bool isfirst) {
	
	//enczData = 
//}




int main(int argc, char **argv) {

		const long dim=64;
	    long sampleDim = dim, factorDim = dim;

		double* zeData = new double[dim*dim];


		if(argc >1){ //인자를 입력했을 경우 해당 경로에 있는 matrix를 불러옴
			string trainfile(argv[1]);
			string testfile=trainfile;
			double* zeData = zDataFromFileFullA(trainfile, factorDim, sampleDim, dim);



		}else{

			for(long j = 0; j < dim; ++j){
				for(long i = 0; i < j ; ++i){
				zeData[j*dim+i]=zeData[i*dim+j]; 	 	
				}

			 	for(long i = j; i < dim; ++i){

				zeData[j*dim+i]= EvaluatorUtils::randomReal(3.0);
				}
			
			} //= random symmetric matrix generation

			
		}


			
		

		//w data generation
			double* wtData = new double[dim];
			 	for(long i = 0; i < dim; ++i){
				wtData[i] = EvaluatorUtils::randomReal(3.0);
				}
		  
		  	double* wData = new double[dim*dim];
			for(long j = 0; j < dim*dim; ++j){
			wData[j] = wtData[j%dim];
			} //horizonal copy generation
		


	long numIter = 1;
	long kdeg = 3;
	double gammaUp = 1;
	double gammaDown = -1;
	bool isInitZero = 0;


	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDim));

	long wBits = 30; 
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
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
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);  //masking matrix - 1st column만 1, 나머지 0 생성
	ZZX poly2 = generateAuxPoly2(slots, batch, pBits, scheme); //masking matrix - 1st row만 1, 나머지 0 생성

	timeutils.stop("Polynomial generation");



	Ciphertext encZData ;
	Ciphertext encWData ;


	timeutils.start("Encrypting Data...");
	encZData = scheme.encrypt(zeData, slots, wBits, logQ);
	encWData = scheme.encrypt(wData, slots, wBits, logQ);

	timeutils.stop("Data encryption");

	timeutils.start("Precomputing");
	Ciphertext AbV= GenAbVertical(encZData, scheme, poly2, bBits, wBits, pBits, slots) ; //각 column이 Ab인 행렬 생성
	Ciphertext AbH= GenAbHorzon(encZData, scheme, poly, bBits, wBits, pBits, slots) ; //각 Row가 Ab인 행렬 생성
	Ciphertext AtA= GenAtA(encZData, scheme, poly, poly2, bBits, wBits, pBits, batch, slots) ; //각 AtA 행렬 생성

	timeutils.stop("Precomputing Done");

		cout << " !!! START ITERATION !!! " << endl;
		
		timeutils.start("mult iter");






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