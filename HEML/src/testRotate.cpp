#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"

#include <cmath>
#include "TrainSVM.h"
#include "CipherSVM.h"

using namespace std;
using namespace NTL;


/* 
 * Test rotate and add! why rotating makes strange result?
 */

int main(int argc, char **argv){

        SetNumThreads(8);

        const long slots = 8;
        double* Avec = new double[slots];
        double* Bvec = new double[slots];


        for(int i =0;i<slots;++i){
                Avec[i] = i+1.0;
                Bvec[i] = -i+8.0;
                cout<<Avec[i]<<", ";
        }
        cout<<endl;

        long wBits = 30;
        long pBits = 20;
        long lBits = 5, aBits = 3;
        long logQ = 400;
        long logN = 14;
        
        Context context(logN,logQ);
        SecretKey secretKey(logN);
        Scheme scheme(secretKey, context);

        scheme.addLeftRotKeys(secretKey);
        scheme.addRightRotKeys(secretKey);
        
        Ciphertext encAvec = scheme.encrypt(Avec,slots,wBits,logQ);
        Ciphertext encBvec = scheme.encrypt(Bvec,slots,wBits, logQ);

        Ciphertext tmp = scheme.rightRotate(encAvec,3);
        
        
