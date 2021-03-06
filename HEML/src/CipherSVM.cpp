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

Ciphertext CipherSVM::GenAtA(Ciphertext encZData, ZZX& poly, ZZX& poly2, long bBits, long wBits, long pBits, long batch, long slots) {
	
	Ciphertext AtA = scheme.multByPoly(encZData, poly2, pBits);
		for (long l = 0; l < bBits; ++l) {
			Ciphertext tmp = scheme.rightRotateByPo2(AtA, l+bBits);//parameter check
			scheme.addAndEqual(AtA, tmp);
		}
                //cout<<"test mult Poly and addAndEqual"<<endl;
                //printDecCiphtxt(scheme.add(AtA,encZData));

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
Ciphertext CipherSVM::GenEncAtA(Ciphertext encZData, ZZX& poly, ZZX& poly2, long bBits, long wBits, long pBits, long batch, long slots,long logQ){
	Ciphertext AtA,tmp,tmp2;
        /*
        complex<double>* zeros = new complex<double>[slots];
        for(long i=0;i<slots;++i){
            zeros[i] =0;
        }
        AtA = scheme.encrypt(zeros,slots,wBits,logQ);
        cout<<"initial AtA"<<endl;
        printDecCiphtxt(AtA);
        //printDecCiphtxt(scheme.add(AtA,AtA));
        //
        */
        //cout<<"print initial AtA.logq: "<<AtA.logq<<endl;

	for(long i=0;i<batch;++i){
		tmp = scheme.leftRotate(encZData,i);
                //cout<<"print initial tmp"<<endl;
                //printDecCiphtxt(tmp);
                //printDecCiphtxt(encZData);
		scheme.multByPolyAndEqual(tmp,poly,pBits);//[a1,0,0;a2,0,0;a3,0,...]
                //scheme.reScaleByAndEqual(tmp,pBits); 
                //cout<<"print tmp"<<endl;
                //printDecCiphtxt(tmp);
                /*
                scheme.modDownToAndEqual(AtA,tmp.logq);
                cout<<"print AtA after add: "<<AtA.logq<<endl;
                scheme.addAndEqual(AtA,tmp);
                printDecCiphtxt(AtA);
                */
		for(long l=0;l<bBits;++l){
			Ciphertext rot = scheme.rightRotateByPo2(tmp,l);
			scheme.addAndEqual(tmp,rot);
	        }
                //cout<<"print tmp after padding"<<endl;
                //printDecCiphtxt(tmp);
                //cout<<"print tmp.logq before Product: "<<tmp.logq<<endl;
		tmp2 = encVerticalVecProduct(encZData,tmp,poly2,bBits,wBits,pBits);
                //cout<<"print tmp2.logq: "<<tmp2.logq<<endl;
                //printDecCiphtxt(tmp2);
		scheme.multByPolyAndEqual(tmp2,poly2,pBits);
                //scheme.reScaleByAndEqual(tmp2,pBits);

                //cout<<"print tmp.logq: "<<tmp.logq<<endl;
		tmp = scheme.rightRotate(tmp2,batch*i);
		//cout<<"print tmp in GenEncAtA"<<endl;
		//printDecCiphtxt(tmp);
                //cout<<"print tmp.q after rotation: "<<tmp.logq<<endl;
                
                //scheme.modDownToAndEqual(AtA,tmp2.logq);
                //cout<<"print AtA.q after modDown: "<<AtA.logq<<endl;
                if(i==0){
                        AtA = tmp;
                        cout<<"0 print"<<endl;
                }else{
                        scheme.addAndEqual(AtA,tmp);
                }
                //printDecCiphtxt(AtA);
		//scheme.addAndEqual(AtA,tmp2);
		//cout<<"print resulting AtA: "<<endl;
                //printDecCiphtxt(AtA);
	}
	return AtA;
}


Ciphertext CipherSVM::GenAbHorzon(Ciphertext encZData,  ZZX& poly, long bBits, long wBits, long pBits, long slots) {
//reesults is [v1,0,...;v2,0,...;...]	
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
        //scheme.reScaleByAndEqual(encIPvec,pBits);
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, l);
			scheme.addAndEqual(encIPvec, rot);
		}


	//scheme.reScaleByAndEqual(encIPvec, wBits);
	scheme.multByPolyAndEqual(encIPvec, poly, pBits);
        //scheme.reScaleByAndEqual(encIPvec,pBits);
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l);
		scheme.addAndEqual(encIPvec, tmp);
	}
	return encIPvec;
	//the result is [v1,v1,...,v1;v2,v2,...,v2;...]
}

Ciphertext CipherSVM::GenAbVertical(Ciphertext encZData,  ZZX& poly, long bBits, long wBits, long pBits, long slots) {
	//results is [v1,v2,v3,...;0,0,...]	
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

	//scheme.reScaleByAndEqual(encIPvec, wBits);
	
		scheme.multByPolyAndEqual(encIPvec, poly, pBits); //Vert poly
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l+bBits);//parameter check
		scheme.addAndEqual(encIPvec, tmp);
	}
	return encIPvec;
}


//GenAbHor

//여기에서의 matrix * vector 연산은 반드시 matrix가 symmetric이여야지 reasonable! maybe...
Ciphertext CipherSVM::encHorizonVecProduct(Ciphertext encZData, Ciphertext encWData, ZZX& poly, long bBits, long wBits, long pBits) {
	Ciphertext encIPvec;
	//original encWData padding should be [w1,w2,...;w1,w2,...]
                cout<<"encZData.logq, encWData.logq: "<<encZData.logq<<", "<<encWData.logq<<endl;
		encIPvec = scheme.modDownTo(encZData, encWData.logq);
		scheme.multAndEqual(encIPvec, encWData); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, l);
			scheme.addAndEqual(encIPvec, rot);
		}
	//the result is c = [c1,x,x,...;c2,x,x,...;c3,x,x,...;...]
	

	scheme.reScaleByAndEqual(encIPvec, wBits);
        /*
	cout<<"before mult poly"<<endl;
        scheme.modDownToAndEqual(encZData,encIPvec.logq);
        printDecCiphtxt(scheme.add(encZData,encIPvec));
	*/
        scheme.multByPolyAndEqual(encIPvec, poly, pBits);
	//scheme.reScaleByAndEqual(encIPvec,pBits);
        for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l);
		scheme.addAndEqual(encIPvec, tmp);
	}
        //cout<<"encZData-encIPvec"<<endl;
        //scheme.modDownToAndEqual(encZData,encIPvec.logq);
        //printDecCiphtxt(scheme.sub(encZData,encIPvec));
	return encIPvec;
	//the result is represented as [v1,v1,v1,..;v2,v2,v2,...]
}

//poly2 V를 넣어야 한다.(poly: 1,1,1,...,1;0,0,0,...;)
Ciphertext CipherSVM::encVerticalVecProduct(Ciphertext encZData, Ciphertext encWData, ZZX& poly,  long bBits, long wBits, long pBits) {
	Ciphertext encIPvec;
	//encWData padding may be [w1,w1,w1,...;w2,w2,w2,...]
		encIPvec = scheme.modDownTo(encZData, encWData.logq);
		scheme.multAndEqual(encIPvec, encWData); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, bBits+l);//paratmeter check
			scheme.addAndEqual(encIPvec, rot);
		}
	
	//the result is c = [c1,c2,c3,...;x,x,x,...;...]

	scheme.reScaleByAndEqual(encIPvec, wBits);
	
	scheme.multByPolyAndEqual(encIPvec, poly, pBits); //Vert poly
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l+bBits);//parameter check
		scheme.addAndEqual(encIPvec, tmp);
	}
	return encIPvec
	;
}
Ciphertext CipherSVM::encVerticalVecProduct2(Ciphertext encZData, Ciphertext encWData, ZZX& poly,  long bBits, long wBits, long pBits) {
	Ciphertext encIPvec;
	//encWData padding may be [w1,w1,w1,...;w2,w2,w2,...]
        // Add print PART!!!
        //
		encIPvec = scheme.modDownTo(encZData, encWData.logq);
		scheme.multAndEqual(encIPvec, encWData); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec, bBits+l);//paratmeter check
			scheme.addAndEqual(encIPvec, rot);
		}
	
	//the result is c = [c1,c2,c3,...;x,x,x,...;...]

	scheme.reScaleByAndEqual(encIPvec, wBits);
        cout<<"print intermedate VerticaVecProduct logp,logq:"<<encIPvec.logp<<", "<<encIPvec.logq<<endl;

	printDecCiphtxt(encIPvec);
	scheme.multByPolyAndEqual(encIPvec, poly, pBits); //Vert poly
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIPvec, l+bBits);//parameter check
		scheme.addAndEqual(encIPvec, tmp);
	}
	return encIPvec
	;
}
ZZX CipherSVM::generateAuxPoly(long slots, long batch, long pBits) {
	//poly is [1,0,0,...;1,0,0,...;...]
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}
ZZX CipherSVM::generateAuxPoly2(long slots, long batch, long pBits) {
	//poly2 is [1,1,1,...,1;0,0,...0;0,0,...]
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < batch; j ++) {   //parameter check
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}
ZZX CipherSVM::generateAuxPolyOne(long slots, long batch, long pBits){
        complex<double>* pvals = new complex<double>[slots];
        for(long j=0;j<slots;++j){
                pvals[j].real(1.0);
        }
        ZZX msg = scheme.context.encode(pvals,slots,pBits);
        delete[] pvals;
        return msg;
}
ZZX CipherSVM::generateAuxPolyConst(double cnst, long slots, long pBits){
    complex<double>* pvals = new complex<double>[slots];
    for(long j=0;j<slots;++j){
        pvals[j].real(cnst);
    }
    ZZX msg = scheme.context.encode(pvals,slots,pBits);
    delete[] pvals;
    return msg;
}

///// Shoudl add the function alternating the direction of padding.

void CipherSVM::encLGDstep(Ciphertext& encWData, Ciphertext& encGrad, double lr,long pBits){
	//NTL_EXEC_RANGE(cnum, first, last);
	//for (long i = first; i < last; ++i) {
        cout<<"initial encGrad.logp: "<<encGrad.logp<<endl;
        printDecCiphtxt(encGrad);
	scheme.multByConstAndEqual(encGrad,lr,pBits);
	cout<<"Print lr*grad (encGrad.logp):"<<encGrad.logp<<endl;
        cout<<"encGrad.logq:"<<encGrad.logq<<endl;
	scheme.reScaleByAndEqual(encGrad,pBits);
	CipherSVM::printDecCiphtxt(encGrad);
	scheme.modDownToAndEqual(encWData, encGrad.logq);
	//should multiply learning rate to gradient!
	scheme.subAndEqual(encWData, encGrad);
	cout<<"Print W-lr*grad"<<endl;
	CipherSVM::printDecCiphtxt(encWData);
	//}
	//NTL_EXEC_RANGE_END;
}
void CipherSVM::encLGDiteration(Ciphertext& encAtAData, Ciphertext& encAbV, Ciphertext& encAbH, Ciphertext& encWData, ZZX& poly,ZZX& poly2,ZZX& polyOne, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits,long slots) {
        
 	//Ciphertext* encGrad = new Ciphertext[cnum];
	//Ciphertext encGrad 
        //generate lr polynomial
        //ZZX polylr = generateAuxPolyConst(gamma, slots, pBits);
	//= new Ciphertext[cnum];	
	cout<<"Initial W (encWData.logq): "<<encWData.logq<<endl;
        cout<<"Initial encWData.lop: "<<encWData.logp<<endl;
	CipherSVM::printDecCiphtxt(encWData);
	//initial Encw = [w1,w2,...;w1,w2,...;...]
        //The important part is pre-computation!
        if(encWData.logq>encAtAData.logq){
            scheme.modDownToAndEqual(encWData,encAtAData.logq);
            cout<<"encWData.logq: "<<encWData.logq<<endl;
        }else{
                scheme.modDownToAndEqual(encAtAData,encWData.logq);
        }
            
        Ciphertext encIP = encHorizonVecProduct(encAtAData, encWData,  poly, bBits, wBits, pBits) ;
	//result is gw = [g1,g1,....,g1;g2,g2,....,g2;....]
        cout<<"Print encAbH.logq: "<<encAbH.logq<<endl;
	cout<<"Print AtA*W (encIP.logq): "<<encIP.logq<<endl;
	CipherSVM::printDecCiphtxt(encIP);
	// 위에 거 결과로 encAtAData와 encWData가 rescale 되지는 않는다.
	scheme.modDownToAndEqual(encAbH,encIP.logq);
        scheme.modDownToAndEqual(encAtAData,encIP.logq);
        cout<<"print encAbH.logp: "<<encAbH.logp<<endl;
        cout<<"print encIP.logp: "<<encIP.logp<<endl;
        cout<<"print encAtAData.logp: "<<encAtAData.logp<<endl;
        //printDecCiphtxt(encAbH);
        scheme.multByPolyAndEqual(encAbH,polyOne,pBits);
        scheme.multByPolyAndEqual(encAbH,polyOne,pBits);
        scheme.reScaleByAndEqual(encAbH,pBits);
        scheme.reScaleByAndEqual(encIP,pBits);
        cout<<"print encAbH.logp: "<<encAbH.logp<<endl;
	scheme.subAndEqual(encIP,encAbH);
	cout<<"Print AtA*W-Atb"<<endl;
	CipherSVM::printDecCiphtxt(encIP);
        //scheme.multByPolyAndEqual(encIP,polyOne,pBits);
	encWData = encHorizonVecProduct(encIMat,encWData,poly,bBits,wBits,pBits);
	cout<<"repadded W (encWData.logp,logq): "<<encWData.logp<<", "<<encWData.logq<<endl;
        cout<<"encIP.logp: "<<encIP.logp<<", encIP.logp: "<<encIP.logq<<endl;
        scheme.multByPolyAndEqual(encWData,polyOne,pBits);
        scheme.multByPolyAndEqual(encWData,polyOne,pBits);
	CipherSVM::printDecCiphtxt(encWData);
	//should repadding the weight vector
	
        encLGDstep(encWData, encIP,gamma,pBits); 

	//scheme.reScaleByAndEqual(encAtAData,wBits);
	//encHorizon...에서 먼저 modDown을 해 주고 시작
	encIP = encVerticalVecProduct2(encAtAData, encWData,  poly2, bBits, wBits, pBits) ;
        cout<<"print 2nd time encIP"<<endl;
        printDecCiphtxt(encIP);
	scheme.modDownToAndEqual(encAbV,encIP.logq);
	scheme.subAndEqual(encIP,encAbV);
	encWData = encVerticalVecProduct(encIMat,encWData,poly2,bBits,wBits,pBits);
        cout<<"error part test"<<endl;
	encLGDstep(encWData, encIP,gamma,wBits); 
	
	//delete[] encIP;//pointer가 아니어서 error???
}
////////Enc Matrix should be obtained!!!!!!////////
void CipherSVM::encMatrix(Ciphertext& encZData, double** zData, long dim, long slots, long batch, long wBits, long logQ){
	complex<double>* pzData = new complex<double>[slots];
	for (long i=0;i<dim;++i){
		for (long j=0;j<dim;++j){
			pzData[i*batch+j].real(zData[i][j]);
		}
	}
	long rest = batch - dim;
	for (long i = 0;i<dim;i++){
		for(long j=dim;j<batch;j++){
			pzData[i*batch+j] = 0;
		}
	}
	for(long i = dim;i<batch;i++){
		for(long j=0;j<batch;j++){
			pzData[i*batch+j] = 0;
		}
	}

	//print pzData
	cout<<"print padded matrix!"<<endl;
	for(long i = 0;i<slots;i++){
		cout<<pzData[i].real()<<", ";
	}
	cout<<endl;
	encZData = scheme.encrypt(pzData,slots,wBits,logQ);
	delete[] pzData;
}
void CipherSVM::printDecCiphtxt(Ciphertext encData){
	complex<double>* dcw = scheme.decrypt(secretKey,encData);
	cout<<"Print ciphertxt but size 9->16"<<endl;
	for (long i = 0;i<16;++i){
		cout<<dcw[i].real()<<", ";
	}
	cout<<endl;
}
//void CipherSVM::decWData(double* wData, Ciphertext encWData, long wBits){}
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

