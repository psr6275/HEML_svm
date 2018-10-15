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


//#include "CipherGD.h"
//#include "TestGD.h"
#include "CipherSVM.h"

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

/// start to define the functions for CipherSVM

/*
// we need to implemet making kernel matrix, but now, we just bring the full A matrix from csv file
// we just encrypt the full A matrix at once.
void CipherSVM::encZData(Ciphertext* encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits, long logQ) {
	complex<double>* pzData = new complex<double>[slots];
	// it should be implemented to deal with large-scale A matrix, 
	//but now we assume that A is not so big.
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < sampleDim; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(zData[j][batch * i + l]);
			}
		}
		encZData[i] = scheme.encrypt(pzData, slots, wBits, logQ);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < sampleDim; ++j) {
		for (long l = 0; l < rest; ++l) {
			pzData[batch * j + l].real(zData[j][batch * (cnum - 1) + l]);
		}
		for (long l = rest; l < batch; ++l) {
			pzData[batch * j + l] = 0;
		}
	}
	encZData[cnum - 1] = scheme.encrypt(pzData, slots, wBits, logQ);

	delete[] pzData;
}
*/

///
Ciphertext CipherSVM::GenAtA(Ciphertext encZData, Scheme& scheme, ZZX& poly, ZZX& poly2, long bBits, long wBits, long pBits, long batch, long slots) {
	
	Ciphertext AtA = scheme.multByPoly(encZData, poly2, pBits);
		for (long l = 0; l < bBits; ++l) {
			Ciphertext tmp = scheme.rightRotateByPo2(AtA, l+bBits);//parameter check
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


Ciphertext CipherSVM::GenAbHorzon(Ciphertext encZData, Scheme& scheme, ZZX& poly, long bBits, long wBits, long pBits, long slots) {
	
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
	return encIPvec;
}

Ciphertext CipherSVM::GenAbVertical(Ciphertext encZData, Scheme& scheme, ZZX& poly, long bBits, long wBits, long pBits, long slots) {
		
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
	return encIPvec;
}


//GenAbHor


Ciphertext CipherSVM::encHorizonVecProduct(Ciphertext encZData, Ciphertext encWData, Scheme& scheme, ZZX& poly, long bBits, long wBits, long pBits) {
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

//poly V를 넣어야 한다.
Ciphertext CipherSVM::encVerticalVecProduct(Ciphertext encZData, Ciphertext encWData, Scheme& scheme, ZZX& poly,  long bBits, long wBits, long pBits) {
	Ciphertext encIPvec;
		encIPvec = scheme.modDownTo(encZData, encWData.logq);
		scheme.multAndEqual(encIPvec, encWData); // xy * w
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
	return encIPvec;
}

ZZX CipherSVM::generateAuxPoly(long slots, long batch, long pBits) {
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}
ZZX CipherSVM::generateAuxPoly2(long slots, long batch, long pBits, Scheme& scheme) {
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < batch; j ++) {   //parameter check
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}

void CipherSVM::encLGDstep(Ciphertext* encWData, Ciphertext* encGrad){
	//NTL_EXEC_RANGE(cnum, first, last);
	//for (long i = first; i < last; ++i) {
	scheme.modDownToAndEqual(encWData[i], encGrad[i].logq);
	scheme.subAndEqual(encWData[i], encGrad[i]);
	//}
	//NTL_EXEC_RANGE_END;
}
void CipherSVM::encLGDiteration(Ciphertext encAtAData, Ciphertext encAbV, Ciphertext encAbH, Ciphertext* encWData, ZZX& poly,ZZX& poly2, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	//Ciphertext* encGrad = new Ciphertext[cnum];
	//Ciphertext encGrad 
	//= new Ciphertext[cnum];	
	Ciphertext encIP = encVerticalVecProduct(encAtAData, encWData, scheme, poly2, bBits, wBits, pBits) ;
	// 위에 거 결과로 encAtAData와 encWData가 rescale 되지는 않는다.
	scheme.modDownToAndEqual(encAbH,encIP.logq);
	scheme.subAndEqual(encIP,encAbH);
	encLGDstep(encWData, encIP, cnum); 

	//scheme.reScaleByAndEqual(encAtAData,wBits);
	//encHorizon...에서 먼저 modDown을 해 주고 시작
	encIP = encHorizonVecProduct(encAtAData, encWData, scheme, poly, bBits, wBits, pBits) ;
	scheme.modDownToAndEqual(encAbV,encIP.logq);
	scheme.subAndEqual(encIP,encAbV);
	encLGDstep(encWData, encGrad, cnum); 
	
	delete[] encIP;
}
void CipherSVM::decWData(double* wData, Ciphertext* encWData, long wBits){}
//////////////////////
//////////////////////
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


Ciphertext makeMatrixA(Ciphertext enczData, Scheme& scheme, long bBits, ZZX& poly, long pBits, long wBits){
		Ciphertext Row;
		Row = scheme.rightRotate(enczData, 1);
		//첫 슬롯에 1 더하기
		for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(Row, l+bBits);//parameter check
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

			const long dim=4;
	/*        long sampleDim = dim, factorDim = dim;

	        	bool isYfirst = atoi(argv[2]);
double** zData;


	if(argc >1){
			long sampleDim = 0, factorDim = 0;
			string trainfile(argv[1]);
			string testfile=trainfile;
			double** zData = zDataFromFile(trainfile, factorDim, sampleDim, isYfirst);

		}else{*/
	        long sampleDim = dim, factorDim = dim;
			double** zData = new double*[dim];
			for(long j = 0; j < dim; ++j){
			double* zj = new double[dim];
				for(long i = 0; i < j ; ++i){
				zj[i]=zData[i][j]; 	
				}

			 	for(long i = j; i < dim; ++i){

				zj[i] = EvaluatorUtils::randomReal(3.0);
				}
			zData[j] = zj;
			} //= random  128 by 128 matrix generation
	//	}


			double* zeData = new double[dim*dim];
			for(long j = 0; j < dim; ++j){
				for(long i = 0; i < dim ; ++i){
				zeData[i*dim+j]=zData[i][j]; 	
				}
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
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	ZZX poly2 = generateAuxPoly2(slots, batch, pBits, scheme);

	timeutils.stop("Polynomial generation");



	Ciphertext encZData ;
	Ciphertext encWData ;


	timeutils.start("Encrypting Data...");
	encZData = scheme.encrypt(zeData, slots, wBits, logQ);
	encWData = scheme.encrypt(wData, slots, wBits, logQ);

	timeutils.stop("Data encryption");

	timeutils.start("Precomputing");
		Ciphertext AbV= GenAbVertical(encZData, scheme, poly2, bBits, wBits, pBits, slots) ; 

	Ciphertext AbH= GenAbHorzon(encZData, scheme, poly, bBits, wBits, pBits, slots) ;
	Ciphertext AtA= GenAtA(encZData, scheme, poly, poly2, bBits, wBits, pBits, batch, slots) ;

	timeutils.stop("Precomputing Done");

		cout << " !!! START ITERATION !!! " << endl;
		
		timeutils.start("mult iter");






//////////one-shot two step

 		Ciphertext encIP = encHorizonVecProduct(encZData, encWData, scheme, poly,  bBits, wBits, pBits) ;

		encWData = encVerticalVecProduct(encZData, encIP, scheme, poly2, bBits, wBits, pBits) ;
			scheme.addAndEqual(encWData, AbH);

		complex<double>* msgg = scheme.decrypt(secretKey, encWData);



//////////////////////

		////Vertical Only
	/*	Ciphertext encIP = encVerticalVecProduct(AtA, encWData, scheme, poly2, bBits, wBits, pBits) ;
		scheme.addAndEqual(encIP, abH);


		complex<double>* msgg = scheme.decrypt(secretKey, encIP);


	
*/


		timeutils.stop("mult iter end");

	


	return 0;
}