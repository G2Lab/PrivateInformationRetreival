/*
Auxiliary functions for integer encoding
*/

#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <helib/helib.h>
#include <helib/Ctxt.h>
#include <helib/polyEval.h>

using namespace std;
using namespace helib;
using namespace NTL;

inline int intlog(unsigned long base, unsigned long input)
{
	int res = max(static_cast<int>(floor(log2(input)/log2(base))),0);
	if(power_long(base, res+1) == input)
		return res+1;
	return res;
}

void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots);

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
void simplePolyEval(Ctxt& ret, const NTL::ZZX& poly, DynamicCtxtPowers& babyStep);

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
void PatersonStockmeyer(Ctxt& ret, const NTL::ZZX& poly, long k, long t, long delta, DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
void degPowerOfTwo(Ctxt& ret, const NTL::ZZX& poly, long k,
        DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

void recursivePolyEval(Ctxt& ret, const NTL::ZZX& poly, long k,
      DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);


// From Geeks for Geeks
// Function for extended Euclidean Algorithm
int gcdExtended(int a, int b, int* x, int* y);
 
// Function to find modulo inverse of a
int modInverse(int A, int M);

int get_inverse(int nom, int dom, int p);

template<typename T, typename Allocator> void print_vector(const vector<T, Allocator>& vect, int num_entries = 10);