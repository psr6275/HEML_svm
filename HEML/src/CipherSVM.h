#ifndef HEML_CIPHERSVM_H_
#define HEML_CIPHERSVM_H_

//#include "NTL/ZZX.h" if we include this file to Library package, we don't have to do this. 
#include <Scheme.h>
#include <SecretKey.h>
#include <complex>
#include "CipherGD.h"

using namespace std;
using namespace NTL;

class CipherSVM {
public:
        Scheme& scheme;
        SecretKey& secretKey;
        Ciphertext encIMat;

        CipherSVM(Scheme& scheme, SecretKey& secretKey, Ciphertext cipherI) : scheme(scheme), secretKey(secretKey), encIMat(cipherI) {}
        //CipherSVM(Scheme& scheme, SecretKey& secretKey);

        void encZData(Ciphertext encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits, long logQ);
        //the first row is 1
        ZZX generateAuxPoly(long slots, long batch, long pBits);
        //the first column is 1
        ZZX generateAuxPoly2(long slots, long batch, long pBits); 
        //generate all 1 polynomial
        ZZX generateAuxPolyOne(long slots, long batch, long pBits);
        ZZX generateAuxPolyConst(double cnst, long slots, long pBits);

        //prepare the basic components
        Ciphertext GenAtA(Ciphertext encZData,  ZZX& poly, ZZX& poly2, long bBits, long wBits, long pBits, long batch, long slots);
        Ciphertext GenAbHorzon(Ciphertext encZData,  ZZX& poly, long bBits, long wBits, long pBits, long slots);
        Ciphertext GenAbVertical(Ciphertext encZData,  ZZX& poly, long bBits, long wBits, long pBits, long slots);
        Ciphertext GenEncAtA(Ciphertext encZData,ZZX& poly, ZZX& poly2, long bBits, long wBits, long pBits,long batch, long slots,long logQ);
        
        //operations for GD step
        Ciphertext encHorizonVecProduct(Ciphertext encZData, Ciphertext encWData,  ZZX& poly, long bBits, long wBits, long pBits); 
        Ciphertext encVerticalVecProduct(Ciphertext encZData, Ciphertext encWData,  ZZX& poly,  long bBits, long wBits, long pBits); 
        Ciphertext encVerticalVecProduct2(Ciphertext encZData, Ciphertext encWData,  ZZX& poly,  long bBits, long wBits, long pBits); 

        //for GD step and iteration
        void encLGDstep(Ciphertext& encWData, Ciphertext& encGrad, double lr,long pBits); 
        void encLGDiteration(Ciphertext& encAtAData,Ciphertext& encAbV, Ciphertext& encAbH, Ciphertext& encWData, ZZX& poly, ZZX& poly2,ZZX& polyOne, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits,long slots);
        //void decWData(double* wData, Ciphertext encWData, long wBits);
        void encMatrix(Ciphertext& encZData, double** zData, long dim, long slots, long bBits, long wBits, long logQ);
        void printDecCiphtxt(Ciphertext encData);

};
#endif
