//max size for the FixedFactor is 32*10 bits
//max chain length is 16
//max sieve size 2^24 = ~16.7M

#define NHALFCHAINLENGTH ((NCHAINLENGTH+1)>>1)

#define ORIGIN_LENGTH 11

#if (ORIGIN_LENGTH != 11)
#error ORIGIN_LENGTH must be 11, other sizes are not supported yet.
#endif

uint GetWordNum(uint nBitNum) {
	return nBitNum >> 5; //32 bits in a uint
}

uint GetBitMask(uint nBitNum) {
	return 1UL << (nBitNum & 31);
}

//vMultipliers is thread-local
void AddMultiplier(__local uint* vMultipliers, const uint nSolvedMultiplier)
{
	uint broken=0;
    // Eliminate duplicates
#pragma unroll
    for (uint i = 0; i < NHALFCHAINLENGTH; i++)
    {
		uint temp = vMultipliers[WORKSIZE*i];
		bool result = !broken && (temp == 0xFFFFFFFF || temp == nSolvedMultiplier);
		vMultipliers[WORKSIZE*i] = select(temp,nSolvedMultiplier,(uint)result);
		broken = select(broken,1U,(uint)result);
    }
}

uint mod(uint a, uint b)
{
	return a%b;
}
uint mulmod(uint a, uint m)
{
	return (((ulong)a)<<32)%m;
}

uint modulo(uint16 n, uint m)
{
	uint result=0;
	result=n.sa%m;

	result=mulmod(result,m);
	result += mod(n.s9,m);
	result = select(result,result-m,result>=m);
	
	result=mulmod(result,m);
	result += mod(n.s8,m);
	result = select(result,result-m,result>=m);
	
	result=mulmod(result,m);
	result += mod(n.s7,m);
	result = select(result,result-m,result>=m);
	
	result=mulmod(result,m);
	result += mod(n.s6,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s5,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s4,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s3,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s2,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s1,m);
	result = select(result,result-m,result>=m);

	result=mulmod(result,m);
	result += mod(n.s0,m);
	result = select(result,result-m,result>=m);

	return result;
}

uint int_invert(uint a, uint nPrime)
{
    int rem0 = nPrime, rem1 = a;
    int aux0 = 0, aux1 = 1;
    int quotient, inverse=0;

	bool go_on = true;
    while (go_on)
    {
		bool test = (rem1<=1);
		inverse = select(inverse,aux1,(int)test);
		go_on = !test;

        quotient = rem0 / rem1;
		rem0 -= quotient * rem1;
        aux0 -= quotient * aux1;

		test = go_on && (rem0 <= 1);
		inverse = select(inverse,aux0,(int)test);
		go_on = select((int)go_on,false,(int)test);
		
        quotient = rem1 / rem0;
		rem1 -= quotient * rem0;
        aux1 -= quotient * aux0;
    }
	return select((uint)inverse,inverse+nPrime,inverse<0);
}

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void CalculateMultipliers(__global const uint*restrict vPrimes, const uint nPrimes, const uint16 nFixedFactor, const uint nCandidatesWords, __global uint* vfCunninghamMultipliers)
{
	__local uint localmem[4*NHALFCHAINLENGTH*WORKSIZE];
	
#pragma unroll
	for(uint i=0; i<4*NHALFCHAINLENGTH*WORKSIZE; i+=WORKSIZE)
		localmem[i+get_local_id(0)] = 0xFFFFFFFF;
	
	__local uint* localCunninghamMultipliers = localmem+get_local_id(0);

	uint nPrimeSeq = get_global_id(0)+1;
    {
        uint nPrime = vPrimes[nPrimeSeq];
        uint nFixedFactorMod = modulo(nFixedFactor, nPrime);
        if (nFixedFactorMod != 0)
		{
			uint nFixedInverse;
			nFixedInverse = int_invert(nFixedFactorMod, nPrime);

			// Weave the sieve for the prime
#pragma unroll
			for (uint nBiTwinSeq = 0; nBiTwinSeq < 2 * NCHAINLENGTH; nBiTwinSeq+=2)
			{
				uint arrnum;
				uint nSolvedMultiplier;
				
				arrnum = select(1,3,nBiTwinSeq+1 >= NCHAINLENGTH);
				nSolvedMultiplier = select(nPrime-nFixedInverse,0U,nFixedInverse==0);
				AddMultiplier(localCunninghamMultipliers+WORKSIZE*arrnum*NHALFCHAINLENGTH, nSolvedMultiplier);

				arrnum = select(0,2,nBiTwinSeq >= NCHAINLENGTH);
				nSolvedMultiplier = nFixedInverse;
				AddMultiplier(localCunninghamMultipliers+WORKSIZE*arrnum*NHALFCHAINLENGTH, nSolvedMultiplier);
				
				//modular divide by 2
				nFixedInverse = select(nFixedInverse>>1,(nFixedInverse+nPrime)>>1,nFixedInverse&1);
			}
			
			//vCunninghamMultipliers+i*ySIZE*xSIZE+nPrimeSeq*xSIZE
			
#pragma unroll
			for(uint h=0; h<4; ++h)
			{
				uint global_base = h*nPrimes*NHALFCHAINLENGTH+nPrimeSeq*NHALFCHAINLENGTH;
				uint local_base = WORKSIZE*h*NHALFCHAINLENGTH;
#pragma unroll
				for(uint i=0; i<NHALFCHAINLENGTH; ++i)
					vfCunninghamMultipliers[global_base+i] = localCunninghamMultipliers[local_base+WORKSIZE*i];
			}
		}
    }
}

#define SIEVE_PER_GROUP (LOCAL_MEM_USED*8)

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void Sieve(__global uint* vfCompositeCunningham, __global const uint* vfCunninghamMultipliers, __global const uint* vPrimes, const uint nPrimes, const uint nCandidatesWords)
{
	__local uint sievepart[SIEVE_PER_GROUP/32];
	
#pragma unroll
	for(uint g=0; g<4; ++g)
	{
#pragma unroll
		for(uint i=0; i<SIEVE_PER_GROUP/32; i+=WORKSIZE)
			sievepart[i+get_local_id(0)] = 0;
		
		uint sievepartition_start=get_group_id(0)*SIEVE_PER_GROUP;
		uint sievepartition_end=get_group_id(0)*SIEVE_PER_GROUP+SIEVE_PER_GROUP;
		sievepartition_end = min(sievepartition_end,(uint)NSIEVESIZE);

		for(uint h=get_local_id(0); h<nPrimes-1; h+=WORKSIZE)
		{
			uint global_base = g*nPrimes*NHALFCHAINLENGTH+(h+1)*NHALFCHAINLENGTH;
			uint nPrime = vPrimes[h+1];
			const uint nRotateBits = nPrime & 31;
			uint nVariableMultiplierStart = sievepartition_start-sievepartition_start%nPrime;
			
			uint broken_mask=0xFFFFFFFF;
#pragma unroll
			for (uint i = 0; i < NHALFCHAINLENGTH; i++)
			{
				uint array_value = vfCunninghamMultipliers[global_base+i];
				broken_mask = select(broken_mask,0U,array_value==0xFFFFFFFF);
				uint nVariableMultiplier = nVariableMultiplierStart+array_value;
				nVariableMultiplier = select(nVariableMultiplier,nVariableMultiplier+nPrime,nVariableMultiplier<sievepartition_start);
				uint lBitMask = GetBitMask(nVariableMultiplier);
				for (; broken_mask!=0 && nVariableMultiplier < sievepartition_end;)
				{
					uint inner_mask=0xFFFFFFFFU;
					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
					inner_mask = select(0U,0xFFFFFFFFU,nVariableMultiplier<sievepartition_end);
					
					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
					inner_mask = select(0U,0xFFFFFFFFU,nVariableMultiplier<sievepartition_end);

					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
					inner_mask = select(0U,0xFFFFFFFFU,nVariableMultiplier<sievepartition_end);

					atomic_or(sievepart+GetWordNum(nVariableMultiplier-sievepartition_start),lBitMask&broken_mask&inner_mask);
					lBitMask = rotate(lBitMask,nRotateBits);
					nVariableMultiplier += nPrime;
				}
			}
		}
#pragma unroll
		for(uint i=0; i<(sievepartition_end-sievepartition_start)/32; i+=WORKSIZE)
			vfCompositeCunningham[g*nCandidatesWords+sievepartition_start/32+i+get_local_id(0)] = sievepart[i+get_local_id(0)];
	}
}

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void Combine(__global const uint* vfCompositeCunningham, const uint nCandidatesWords, __global uint*restrict counter, __global uint*restrict CandidateList)
{
	uint i=get_global_id(0);
	uint n0=vfCompositeCunningham[0*nCandidatesWords+i];
	uint n1=vfCompositeCunningham[1*nCandidatesWords+i];
	uint n2=vfCompositeCunningham[2*nCandidatesWords+i];
	uint n3=vfCompositeCunningham[3*nCandidatesWords+i];
	uint bits = ~((n0|n2) & (n1|n3) & (n0|n1));

	//get rid of the candidate 0. the ALU is way underused in this kernel, so doing this doesn't slow anything down. i love tricks like this :)
	bits &= select(~0,~1,i==0);

#pragma unroll
	for(uint j=0; j<32; ++j)
	{
		if (bits&(1<<j))
		{
			uint counter_value = atomic_inc(counter);
			CandidateList[counter_value] = i*32+j;
		}
	}
}

//Fermat test code:
#define LoopDown(x) x(f); x(e); x(d); x(c); x(b); x(a); x(9); x(8); x(7); x(6); x(5); x(4); x(3); x(2); x(1); x(0); 
#define LoopUp(x) x(0); x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9); x(a); x(b); x(c); x(d); x(e); x(f); 

#define CompareIter(n) \
	highertmp |= select(0U,1U<<0x ## n,a.s ## n>b.s ## n);\
	lowertmp |= select(0U,1U<<0x ## n,a.s ## n<b.s ## n);


bool Compare(uint16 a, uint16 b, const bool lower, const bool equal, const bool higher)
{
	uint lowertmp=0,highertmp=0;
	
	LoopDown(CompareIter);
	bool result = select(select((uint)lower,(uint)higher,highertmp>lowertmp),(uint)equal,lowertmp==highertmp);
	
	return result;
}

bool a_LT_b(uint16 a, uint16 b)
{
	return Compare(a,b,true,false,false);
}
bool a_EQ_b(uint16 a, uint16 b)
{
	return Compare(a,b,false,true,false);
}
bool a_GT_b(uint16 a, uint16 b)
{
	return Compare(a,b,false,false,true);
}

uint GetLog2(uint v)
{
	uint r; // result of log2(v) will go here
	uint shift;
	r = select(0,1<<4,v > 0xFFFF); v >>= r;
	shift = select(0,1<<3,v > 0xFF); v >>= shift; r |= shift;
	shift = select(0,1<<2,v > 0xF ); v >>= shift; r |= shift;
	shift = select(0,1<<1,v > 0x3 ); v >>= shift; r |= shift;
	r |= (v >> 1);
	return r;
}

#define GetBitLengthIter(n) \
	bits = select(bits,a.s ## n,(uint)(a.s ## n != 0 && !found));\
	length = select(length-32U,length,(uint)found);\
	found = select((uint)found,(uint)true,(uint)(a.s ## n != 0));

uint GetBitLength(uint16 a)
{
	uint length=512;
	uint bits=0;
	
	bool found=false;
	
	LoopDown(GetBitLengthIter);
	return length+GetLog2(bits)+1;
}

uint16 a_MINUS1_fast(uint16 a)
{
	a.s0 -= 1;
	return a;
}

uint GetIndex(const uint16 a, uint i)
{
	uint8 b = select(a.s01234567,a.s89abcdef,(uint8)(select(0U,0xFFFFFFFFU,i&8)));
	uint4 c = select(b.s0123    ,b.s4567    ,(uint4)(select(0U,0xFFFFFFFFU,i&4)));
	uint2 d = select(c.s01      ,c.s23      ,(uint2)(select(0U,0xFFFFFFFFU,i&2)));
	uint  e = select(d.s0       ,d.s1       ,       (select(0U,0xFFFFFFFFU,i&1)));
	return e;
}

uint16 SetIndex(const uint a, uint i)
{
	uint16 ret;
	i&=15;
	ret.s0 = select(0U,a,i==0x0);
	ret.s1 = select(0U,a,i==0x1);
	ret.s2 = select(0U,a,i==0x2);
	ret.s3 = select(0U,a,i==0x3);
	ret.s4 = select(0U,a,i==0x4);
	ret.s5 = select(0U,a,i==0x5);
	ret.s6 = select(0U,a,i==0x6);
	ret.s7 = select(0U,a,i==0x7);
	ret.s8 = select(0U,a,i==0x8);
	ret.s9 = select(0U,a,i==0x9);
	ret.sa = select(0U,a,i==0xa);
	ret.sb = select(0U,a,i==0xb);
	ret.sc = select(0U,a,i==0xc);
	ret.sd = select(0U,a,i==0xd);
	ret.se = select(0U,a,i==0xe);
	ret.sf = select(0U,a,i==0xf);
	return ret;
}

uint lessthan(ulong a, ulong b)
{
	bool lt_hi = (a>>32)<(b>>32);
	bool lt_lo = (a&0xFFFFFFFF)<(b&0xFFFFFFFF);
	bool eq_hi = ((a>>32)==(b>>32));
	bool ret = select((uint)lt_hi,(uint)lt_lo,(uint)eq_hi);
	return (uint)ret;
}

#define InnerMul(result,n,minj,maxj) \
	for(int j=minj; j<maxj; ++j)\
	{\
		ulong oldtmp = tmp;\
		tmp += ((ulong)(GetIndex(b,0x0 ## n-j)))*GetIndex(a,j);\
		tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
	}\
	(*result).s ## n =(uint)tmp;\
	tmp = upsample(tmp_higher,tmp>>32);\
	tmp_higher=0;

void Mul(uint16* hi, uint16* lo, uint16 a, uint16 b)
{
	ulong tmp=0;
	
	const int alength=ORIGIN_LENGTH+1;
	const int blength=ORIGIN_LENGTH+1;
	
	const int combined_length = alength+blength;
	uint tmp_higher=0;

#pragma unroll
	InnerMul(lo,0,0,1);
#pragma unroll
	InnerMul(lo,1,0,2);
#pragma unroll
	InnerMul(lo,2,0,3);
#pragma unroll
	InnerMul(lo,3,0,4);
#pragma unroll
	InnerMul(lo,4,0,5);
#pragma unroll
	InnerMul(lo,5,0,6);
#pragma unroll
	InnerMul(lo,6,0,7);
#pragma unroll
	InnerMul(lo,7,0,8);
#pragma unroll
	InnerMul(lo,8,0,9);
#pragma unroll
	InnerMul(lo,9,0,10);
#pragma unroll
	InnerMul(lo,a,0,11);
#pragma unroll
	InnerMul(lo,b,0,12);
#pragma unroll
	InnerMul(lo,c,1,12);
#pragma unroll
	InnerMul(lo,d,2,12);
#pragma unroll
	InnerMul(lo,e,3,12);
#pragma unroll
	InnerMul(lo,f,4,12);

#pragma unroll
	InnerMul(hi,0,5,12);
#pragma unroll
	InnerMul(hi,1,6,12);
#pragma unroll
	InnerMul(hi,2,7,12);
#pragma unroll
	InnerMul(hi,3,8,12);
#pragma unroll
	InnerMul(hi,4,9,12);
#pragma unroll
	InnerMul(hi,5,10,12);
#pragma unroll
	InnerMul(hi,6,11,12);
#pragma unroll
	InnerMul(hi,7,12,12);
}
void Mul_lower(uint16* lo, uint16 a, uint16 b)
{
	ulong tmp=0;
	uint tmp_higher=0;
	
#pragma unroll
	InnerMul(lo,0,0,1);
#pragma unroll
	InnerMul(lo,1,0,2);
#pragma unroll
	InnerMul(lo,2,0,3);
#pragma unroll
	InnerMul(lo,3,0,4);
#pragma unroll
	InnerMul(lo,4,0,5);
#pragma unroll
	InnerMul(lo,5,0,6);
#pragma unroll
	InnerMul(lo,6,0,7);
#pragma unroll
	InnerMul(lo,7,0,8);
#pragma unroll
	InnerMul(lo,8,0,9);
#pragma unroll
	InnerMul(lo,9,0,10);
#pragma unroll
	InnerMul(lo,a,0,11);
#pragma unroll
	InnerMul(lo,b,0,12);
}

#define InnerSquare_even(result,n,minj,maxj) \
	for(int j=minj; 2*j<maxj; ++j)\
	{\
		ulong num = ((ulong)(GetIndex(a,maxj-j)))*GetIndex(a,j);\
		ulong oldtmp = tmp;\
		tmp += num;\
		tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
		oldtmp = tmp;\
		tmp += num;\
		tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
	}\
	{\
		ulong oldtmp = tmp;\
		uint num = GetIndex(a,maxj>>1);\
		tmp += ((ulong)num)*num;\
		tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
	}\
	(*result).s ## n =(uint)tmp;\
	tmp = upsample(tmp_higher,tmp>>32);\
	tmp_higher=0;

#define InnerSquare_odd(result,n,minj,maxj) \
	for(int j=minj; 2*j<maxj; ++j)\
	{\
		ulong num = ((ulong)(GetIndex(a,maxj-j)))*GetIndex(a,j);\
		ulong oldtmp = tmp;\
		tmp += num;\
		tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
		oldtmp = tmp;\
		tmp += num;\
		tmp_higher += select(0U,1U,lessthan(tmp,oldtmp));\
	}\
	(*result).s ## n =(uint)tmp;\
	tmp = upsample(tmp_higher,tmp>>32);\
	tmp_higher=0;

	
//A specialized squaring algorithm should be faster than using a generic multiplication algorithm.
//It just so happens that at least on a HD6990, it isn't. It's actually a bit slower.
//Maybe my implementation is bad? For now, you can toggle using Square() with the define USE_SQUARE_ALGORITHM.
//-mtrlt
void Square(uint16* hi, uint16* lo, uint16 a)
{
	ulong tmp=0;
	uint tmp_higher=0;	
#pragma unroll
	InnerSquare_even(lo,0,0,0);
#pragma unroll
	InnerSquare_odd(lo,1,0,1);
#pragma unroll
	InnerSquare_even(lo,2,0,2);
#pragma unroll
	InnerSquare_odd(lo,3,0,3);
#pragma unroll
	InnerSquare_even(lo,4,0,4);
#pragma unroll
	InnerSquare_odd(lo,5,0,5);
#pragma unroll
	InnerSquare_even(lo,6,0,6);
#pragma unroll
	InnerSquare_odd(lo,7,0,7);
#pragma unroll
	InnerSquare_even(lo,8,0,8);
#pragma unroll
	InnerSquare_odd(lo,9,0,9);
#pragma unroll
	InnerSquare_even(lo,a,0,10);
#pragma unroll
	InnerSquare_odd(lo,b,0,11);
#pragma unroll
	InnerSquare_even(lo,c,1,12);
#pragma unroll
	InnerSquare_odd(lo,d,2,13);
#pragma unroll
	InnerSquare_even(lo,e,3,14);
#pragma unroll
	InnerSquare_odd(lo,f,4,15);

#pragma unroll
	InnerSquare_even(hi,0,5,16);
#pragma unroll
	InnerSquare_odd(hi,1,6,17);
#pragma unroll
	InnerSquare_even(hi,2,7,18);
#pragma unroll
	InnerSquare_odd(hi,3,8,19);
#pragma unroll
	InnerSquare_even(hi,4,9,20);
#pragma unroll
	InnerSquare_odd(hi,5,10,21);
#pragma unroll
	InnerSquare_even(hi,6,11,22);
#pragma unroll
	InnerSquare_odd(hi,7,12,23);

}

uint16 Mul_512b_32b_ORIGIN_LENGTH(const uint16 a, uint b)
{
	ulong overflow=0;
	uint16 ret=0;
#pragma unroll
	for(uint i=0; i<ORIGIN_LENGTH; ++i)
	{
		overflow += ((ulong)(GetIndex(a,i)))*b;
		ret |= SetIndex((uint)overflow,i);
		overflow >>= 32;
	}
	ret |= SetIndex((uint)overflow,ORIGIN_LENGTH);
	return ret;
}

uint16 Mul_512b_32b_ORIGIN_LENGTH_PLUS1(const uint16 a, uint b)
{
	ulong overflow=0;
	uint16 ret=0;
#pragma unroll
	for(uint i=0; i<ORIGIN_LENGTH+1; ++i)
	{
		overflow += ((ulong)(GetIndex(a,i)))*b;
		ret |= SetIndex((uint)overflow,i);
		overflow >>= 32;
	}
	ret |= SetIndex((uint)overflow,ORIGIN_LENGTH+1);
	return ret;
}

bool Mul_512b_32b_compare(uint16 a, uint b, const uint16 compared)
{
	a = Mul_512b_32b_ORIGIN_LENGTH_PLUS1(a,b);
	bool result = !a_LT_b(a,compared);
	return result;
}

uint16 Sub(uint16 a, uint16 b)
{
	uint overflow=0, overflow_old=0;
	uint16 ret=0;
	
#pragma unroll
	for(uint i=0; i<ORIGIN_LENGTH+2; ++i)
	{
		uint a_i = a.s0;
		uint b_i = b.s0;
		overflow_old = overflow;
		overflow = select(0U,1U,a_i<b_i+overflow_old);
		a_i -= b_i+overflow_old;
		ret |= SetIndex(a_i,i);
		a = a.s123456789abcdef0;
		b = b.s123456789abcdef0;
	}
	return ret;
}

//Calculating the Mu value for Barrett Reduction
uint16 Calculate_Mu(uint16 m, uint bits)
{
	uint onebit = bits+2; //the number of the bit that should be 1
	uint16 k = 0;
	k.s8 = 1<<(onebit&31);
	uint16 mu = 0;

	for(int i=(onebit-8*32)>>5; i>=0; --i)
	{
		uint x=0;
		//now k >= m
		//find highest x such that m*x<=k
		for(uint j=0; j<32; ++j)
		{
			uint coef = 0x80000000U>>j;
			bool gt_eq = Mul_512b_32b_compare(m,x|coef,k);
			x = select(x|coef,x,(uint)gt_eq);
		}
		mu |= SetIndex(x,i);
		k = Sub(k,Mul_512b_32b_ORIGIN_LENGTH_PLUS1(m,x));
		k = k.sf0123456789abcde;
	}
	return mu;
}

uint16 SHR_big(uint16 hi, uint16 lo, uint n)
{
	uint n_8 = (n&256);
	lo = select(lo,lo.s89abcdef01234567,(uint16)select(0U,0xFFFFFFFFU,n_8));
	hi = select(hi,hi.s89abcdef01234567,(uint16)select(0U,0xFFFFFFFFU,n_8));
	lo.s89abcdef = select(lo.s89abcdef,hi.s89abcdef,(uint8)select(0U,0xFFFFFFFFU,n_8));

	uint n_4 = (n&128);
	lo = select(lo,lo.s456789abcdef0123,(uint16)select(0U,0xFFFFFFFFU,n_4));
	hi = select(hi,hi.s456789abcdef0123,(uint16)select(0U,0xFFFFFFFFU,n_4));
	lo.scdef = select(lo.scdef,hi.scdef,(uint4)select(0U,0xFFFFFFFFU,n_4));

	uint n_2 = (n&64);
	lo = select(lo,lo.s23456789abcdef01,(uint16)select(0U,0xFFFFFFFFU,n_2));
	hi = select(hi,hi.s23456789abcdef01,(uint16)select(0U,0xFFFFFFFFU,n_2));
	lo.sef = select(lo.sef,hi.sef,(uint2)select(0U,0xFFFFFFFFU,n_2));

	uint n_1 = (n&32);
	lo = select(lo,lo.s123456789abcdef0,(uint16)select(0U,0xFFFFFFFFU,n_1));
	hi = select(hi,hi.s123456789abcdef0,(uint16)select(0U,0xFFFFFFFFU,n_1));
	lo.sf = select(lo.sf,hi.sf,(uint )select(0U,0xFFFFFFFFU,n_1));
		
	n&=31;
	uint n_0 = (n>0);
	
	uint losf_tmp = lo.sf;
	
	lo = select(lo,(lo>>n) | (lo.s123456789abcdef0 << 32-n),(uint16)select(0U,0xFFFFFFFFU,n_0));
	lo.sf = select(lo.sf,(losf_tmp >> n) | (hi.s0 << 32-n),n_0);
	return lo;
}

uint GetFirstBitsMask(const uint bits, const uint wordnum)
{
	const uint word_start_bit = wordnum*32;
	const uint word_end_bit = wordnum*32+31;
	
	uint belowmask = 0;
	uint abovemask = 0xFFFFFFFFU;
	uint betweenmask = 0xFFFFFFFFU>>(32-(bits-word_start_bit));
	
	uint mask = betweenmask;
	mask = select(mask,belowmask,bits <= word_start_bit);
	mask = select(mask,abovemask,bits > word_end_bit);
	return mask;
}

#define GetFirstBitsMasksIter(n) \
	fbmask.s ## n = GetFirstBitsMask(k+1, 0x0 ## n);

uint16 Mod(uint16 hi, uint16 lo, uint16 modulus, uint modulus_length, uint16 mu)
{
	uint k = modulus_length;
	
	uint16 q = SHR_big(hi,lo,k);
	
	uint16 fbmask = 0xFFFFFFFF;
	GetFirstBitsMasksIter(a);
	GetFirstBitsMasksIter(b);
	GetFirstBitsMasksIter(c);
	fbmask.sdef = 0;
	
	uint16 hi2=0,lo2=0;
	Mul(&hi2,&lo2,q,mu);
	q = SHR_big(hi2,lo2,k+2);
	uint16 tmp1 = lo&fbmask;
	
	lo2=0;
	Mul_lower(&lo2,q,modulus);
	uint16 tmp2 = lo2&fbmask;

	uint16 x = 0;
	
	bool altb = a_LT_b(tmp1,tmp2);
	
	tmp2 = select(tmp2,Sub(tmp2,tmp1),(uint16)select(0U,0xFFFFFFFFU,(uint)altb));
	x = SetIndex(select(0,1<<((k+1)&31),(uint)altb),(k+1)>>5);
	
	uint16 selected = select(tmp1,x,(uint16)(select(0U,0xFFFFFFFFU,(uint)altb)));
	x = Sub(selected,tmp2);
	x = select(Sub(x,modulus),x,(uint16)select(0U,0xFFFFFFFFU,(uint)a_LT_b(x,modulus)));
	x = select(Sub(x,modulus),x,(uint16)select(0U,0xFFFFFFFFU,(uint)a_LT_b(x,modulus)));
	return x;
}

uint16 MulMod(uint16 a, uint16 b, uint16 m, uint modulus_length, uint16 mu)
{
	uint16 hi=0,lo=0;
	Mul(&hi, &lo, a,b);
	return Mod(hi,lo,m,modulus_length,mu);
}
uint16 SquareMod(uint16 a, uint16 m, uint modulus_length, uint16 mu)
{
	uint16 hi=0,lo=0;
#ifdef USE_SQUARE_ALGORITHM
	Square(&hi, &lo, a);
#else
	Mul(&hi, &lo, a,a);
#endif
	return Mod(hi,lo,m,modulus_length,mu);
}

uint FermatTest(uint16 modulus)
{
	uint16 result = (uint16)(1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0);
	
	uint16 base = (uint16)(2,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0);
	uint16 exponent = a_MINUS1_fast(modulus);
	
	uint modulus_length = GetBitLength(modulus);
	uint16 mu = Calculate_Mu(modulus,modulus_length*2);

	uint exponent_size = ORIGIN_LENGTH+1;
	
	for(uint i=0; i<exponent_size; ++i)
	{
		uint current_exponent = exponent.s0;
		exponent = exponent.s123456789abcdef0;
		for(uint b=0; b<32; ++b)
		{
			uint16 tmp = MulMod(result,base,modulus,modulus_length,mu);
			result = select(result,tmp,(uint16)(select(0U,0xFFFFFFFFU,(current_exponent>>b)&1)));
			base = SquareMod(base,modulus,modulus_length,mu);
		}
	}
	return result.s0; //good enough for our purposes, results in extremely rare false positives when the Fermat Test remainder is actually of the form X*2^32+1 where X is a positive integer.
}

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void Fermat(__global const uint*restrict FermatNumbers, __global uint*restrict FermatOutput, const uint16 FixedMultiplier)
{
	uint16 a = Mul_512b_32b_ORIGIN_LENGTH(FixedMultiplier,FermatNumbers[get_global_id(0)>>1]);

	//at this point, a has the chain origin. even threads subtract one, odd threads add one. this way primes on both sides of the origin are tested.
	uint added=select(1U,0xFFFFFFFFU,(get_global_id(0)&1) == 0);
	uint compared=select(0U,0xFFFFFFFFU,(get_global_id(0)&1) == 0);
	a.s0 += added;
	a.s1 += select(0U,added,a.s0==compared);
	a.s2 += select(0U,added,a.s1==compared);
	//could add more digits, I suppose. won't really change anything.

	uint ret = FermatTest(a);
	FermatOutput[get_global_id(0)] = ret;
}
