#ifndef HEML_CIPHERSVM_H_
#define HEML_CIPHERSVM_H_

//#include "NTL/ZZX.h" if we include this file to Library package, we don't have to do this. 
#include <Scheme.h>
#include <SecretKey.h>
#include <CipherGD.h>
#include <complex>

using namespace std;
using namespace NTL;

class CipherSVM {
public:
        Scheme& scheme;
        SecretKey& secretKey;

        CipherSVM(Scheme& scheme, SecretKey& secretKey) : scheme(scheme), secretKey(secretKey) {}
        
}
