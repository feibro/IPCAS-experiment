extern "C" {
#include"miracl.h"
}
#include<iostream>
#include <cstdio>
#include <ctime>
#include"big.h"
#include "pairing_1.h"
#define MR_PAIRING_SSP  
#define SECURITY_PARAMETER 80  
Big randbits(int n) { Big z; bigbits(n, z.fn); return z; }
int main() {
    PFC pfc(SECURITY_PARAMETER);  // initialise pairing-friendly curve

    G1 P, Q, R;
    GT T , S ;
    Big s,r;
    char str[64] = { 0 };
    int len = 5000;
    double ave = 0, start = 0, end = 0, all = 0;
    //cout <<  pfc.order() << endl;

    //Pairing Time
    for (int i = 0; i < len; i++)
    {
        pfc.random(P);
        pfc.random(Q);
        start = clock();
        T = pfc.pairing(P,Q);
        end = clock();
        all += end - start;
    }
    ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
    printf("----Pairing Opertion Consuming Time : %.6f ms ----\n", ave);
    
    //Pairing-based Add Test G1
    ave = 0, start = 0, end = 0, all = 0;
    for (int i = 0; i < len; i++)
    {
        pfc.random(P);
        pfc.random(Q);
        start = clock();
        R = P + Q;
        end = clock();
        all += end - start;
    }
    ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
    printf("----Pairing Addition Opertion Consuming Time On G1 : %.6f ms ----\n", ave);

    // Pairing-based Mul Test G1
    ave = 0, start = 0, end = 0, all = 0;
    for (int i = 0; i < len; i++)
    {
        pfc.random(P);
        pfc.random(s);
        start = clock();
        R = pfc.mult(P,s);
        end = clock();
        all += end - start;
    }
    ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
    printf("----Pairing Scalar Multiplication Opertion Consuming Time On G1: %.6f ms ----\n", ave);

    // Pairing-based scalar Pow on GT
    ave = 0, start = 0, end = 0, all = 0;
    for (int i = 0; i < len; i++)
    {
        pfc.random(P);
        pfc.random(Q);
        T = pfc.pairing(P,Q);
        pfc.random(P);
        pfc.random(Q);
        S = pfc.pairing(P, Q);
        start = clock();
        S = S * T;
        end = clock();
        all += end - start;
    }
    ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
    printf("---- Pairing Multiplication Opertion Consuming Time On GT: %.6f ms ----\n", ave);
    
    // map-to-point test
    ave = 0, start = 0, end = 0, all = 0;
    for (int i = 0; i < len; i++)
    {
        pfc.random(P);
        r = randbits(192);
        to_binary(r, 64, str);
        start = clock();
        pfc.hash_and_map(P, str);
        end = clock();
        all += end - start;
    }
    ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
    printf("---- Map To Point Opertion Consuming Time: %.6f ms ----\n", ave);
    
    return 0;
}
