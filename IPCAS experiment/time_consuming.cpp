#define _CRT_SECURE_NO_WARNINGS 
extern "C" {
#include "miracl.h"
}
#include <iostream>
#include <cstdio>
#include <ctime>
#include "big.h"
#include "ecn.h"

//! NIST p192 bits ECC curve prime
char* ecp = (char*)"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF";

//! NIST p192 bits ECC curve parameter b
char* ecb = (char*)"64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1";

//! NIST p192 bits ECC curve parameter q
char* ecq = (char*)"FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831";

//! NIST p192 bits ECC curve point of prime order (x,y)
char* ecx = (char*)"188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
char* ecy = (char*)"07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";

#define HASH_LEN 32
Miracl precision(196,16);

Big H1(const char* string)
{ // Hash a zero-terminated string to a number < modulus
	Big h;
	char s[HASH_LEN];
	int i;
	sha256 sh;
	shs256_init(&sh);

	for (i = 0;; i++)
	{
		if (string[i] == 0) break;
		shs256_process(&sh, string[i]);
	}
	shs256_hash(&sh, s);
	//h = from_binary(HASH_LEN, s);
	return h;
}
int main()
{
	time_t seed;
	Big a, b, p, x, y, q, s, hash, r;
	ECn P, point1, point2, point3;
	double ave = 0, start = 0, end = 0, all = 0;
	int len = 5000;

	time(&seed);
	irand((long)seed);

	miracl* mip = &precision;

	// ECC init
	a = -3;
	mip->IOBASE = 16;
	b = ecb;
	p = ecp;
	q = ecq;
	ecurve(a, b, p, MR_BEST);

	x = ecx;
	y = ecy;
	P = ECn(x, y);

	//ECC Addition Time
	start = end = all = ave = 0;
	for (int i = 1; i < len; i++)
	{
		point1 = rand(q) * P;
		point2 = rand(q) * P;
		start = clock();
		point1 += point2;
		end = clock();
		all += end - start;
	}
	ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
	printf("----ECC Addition Operation Consuming Time : %.6f ms ----\n", ave);

	//ECC Multiplication Time
	start = end = all = ave = 0;
	for (int i = 1; i < len; i++)
	{
		s = rand(q);
		start = clock();
		point3 = s * P;
		end = clock();
		all += end - start;
	}
	ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
	printf("----ECC Multiplication Operation Consuming Time : %.6f ms ----\n", ave);

	//Sha-256 Hash Time
	start = end = all = ave = 0;
	for (int i = 0; i < len; i++) {
		char str[64] = { 0 };
		r = rand(q);
		to_binary(r, 64, str);
		start = clock();
		hash = H1(str);
		end = clock();
		all += end - start;
	}
	ave = all / (double)len / CLOCKS_PER_SEC * 1000.0;
	printf("----Hash Operation Consuming Time : %.6f ms ----\n", ave);
	
	return 0;
}


