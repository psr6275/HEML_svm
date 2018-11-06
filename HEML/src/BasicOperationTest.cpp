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


int main(){
    const long dim=4;
    long m = 3,n=2;
    long fdimBits = (long)ceil(log2(n));
    long sdimBits = (long)ceil(log2(m));
    long sBits = fdimBits+ sdimBits;
    long batch = 1<<fdimBits;//dimensoin을 2의 power로 다시 구한것! 넉넉히 더 효율적인 벡터 메트릭스 곱 위해서
    long slots = 1<<sBits;
    double* A = new double[m*n];
    double* b = new double[n];
    b[0] = 1;
    b[1] = 3;
    double* B = new double[m*n];
    //make a A matrix!
    for(long i = 0; i<m;i++){
        for(long j=0;j<n;j++){
            A[i*m+j] = i*m+j;
            B[i*m+j] = b[j];
            cout<< "idx = "<<i*m+j<<", A[idx] = "<<A[i*m+j]<<", B[idx] = "<<B[i*m+j]<<endl;
        }
    }
    long logN=14;
    long logQ=400;
    long pBits=40;//polynomial encode 와 관련 있는 것으로 왠지 scaling factor와 관련 있는듯!
    //근데 그냥 바로 real-valued 부터 encrypt하면 그런건 어떻게 설정되지?
    long wBits=30;//logp 즉 한번 rescale시 깎이는 것고 관련!

    Context context(logN,logQ);
    SecretKey secretKey(logN);
    Scheme scheme(secretKey,context);

    complex<double>* DataA = new complex<double>[slots];
    complex<double>* DataB = new complex<double>[slots];
    Ciphertext cipherA, cipherB;
    for(long i = 0;i<m;i++){
        for(long j = 1;j<n;j++){
            DataA[batch*i+j].real(A[m*i+j]);
            DataB[batch*i+j].real(b[j]);
        }
    }

    for(long l=0;l<slots;l++){
        cout<<"print DataA: "<<DataA[l]<<" | print DataB: "<<DataB[l]<<endl;
    }
    
    cipherA = scheme.encrypt(DataA,slots,wBits, logQ);
    cipherB = scheme.encrypt(DataB,slots,wBits,logQ)

    complex<double>* pvals = new complex<double>[slots];
    for(long i=0;i<slots;i+=batch){
        pvals[i].real(1.0);
    }
    ZZX msg = scheme.context.encode(pvals,slots,pBits);
    ZZX msg2 = scheme.encode(pvals,slots,pBits);
    cout<<"msg: "<<msg<<endl;
    cout<<"msg2: "<<msg2<<endl;
    //cout<<"msg2: logp "<<msg2.logp<<", logq "<<msg2.logq<<endl; 

    Ciphertext encIP = scheme.mult(cipherA,cipherB);
    
    for(long l = 0;l<fdimBits;l++){
        Ciphertext rot = scheme.leftRotateByPo2(encIP,l);
        scheme.addAndEqual(encIP,rot);    
    }
    cout<<"before rescale encIP.logq: "<<encIP.logq<<endl;
    scheme.reScaleByAndEqual(encIP,wBits);
    cout<<"after rescale encIP.logq: "<<encIP.logq<<endl;

    scheme.multByPolyAndEqual(encIP,msg,pBits);
    //cout<<"msg: logp "<<msg.logp<<", logq "<<msg.logq<<endl;
    cout<<"after polynomial mult encIP.logq: "<<encIP.logq<<endl;

    complex<double>* DResults = scheme.decrypt(secretKey,encIP);
    cout<<"Results: "<<endl;
    for (long i = 0; i<slots;i++){
        cout<< DResults[i] <<endl;
    }
    
    return 0;




    
}